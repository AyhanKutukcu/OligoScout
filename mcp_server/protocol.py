"""Dependency-free Model Context Protocol (MCP) server core for OligoScout.

Implements the JSON-RPC 2.0 message layer of MCP (revisions 2025-06-18,
2025-03-26 and 2024-11-05) with two transports:

* stdio            -- newline-delimited JSON on stdin/stdout (Claude Desktop,
                      Claude Code, Cursor, VS Code, Gemini CLI, LM Studio, ...)
* Streamable HTTP  -- POST /mcp returning application/json (ChatGPT connectors,
                      claude.ai remote connectors, any HTTP MCP client)

The same HTTP server also exposes every tool as a plain REST endpoint
(POST /api/tools/<name>) with an OpenAPI 3.1 document at /openapi.json for
platforms that use function calling instead of MCP (GPT Actions, Gemini
function declarations, generic agents).

Only the Python standard library is used, matching the reference engine.
"""
from __future__ import annotations

import hmac
import json
import os
import secrets
import sys
import threading
import time
import traceback
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

from oauth import parse_form

SUPPORTED_VERSIONS = ("2025-06-18", "2025-03-26", "2024-11-05")
LATEST = SUPPORTED_VERSIONS[0]

PARSE_ERROR, INVALID_REQUEST, METHOD_NOT_FOUND, INVALID_PARAMS, INTERNAL_ERROR = (
    -32700, -32600, -32601, -32602, -32603)


class ToolInputError(ValueError):
    """Raised for invalid tool arguments; surfaced as a tool error, not a protocol error."""


def log(*parts):
    print("[oligoscout-mcp]", *parts, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Minimal JSON-Schema validation (the subset used by our tool schemas)
# ---------------------------------------------------------------------------
_TYPES = {"string": str, "integer": int, "number": (int, float), "boolean": bool,
          "object": dict, "array": list}


def validate(value, schema, path="arguments"):
    t = schema.get("type")
    if t:
        types = t if isinstance(t, list) else [t]
        ok = False
        for name in types:
            if name == "null" and value is None:
                ok = True
            elif name in _TYPES and isinstance(value, _TYPES[name]) and not (
                    name in ("integer", "number") and isinstance(value, bool)):
                ok = not (name == "integer" and isinstance(value, float) and not value.is_integer())
            if ok:
                break
        if not ok:
            raise ToolInputError(f"{path}: beklenen tür {t}, gelen {type(value).__name__}")
    if "enum" in schema and value not in schema["enum"]:
        raise ToolInputError(f"{path}: {value!r} geçersiz; izinli: {schema['enum']}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            raise ToolInputError(f"{path}: en az {schema['minimum']} olmalı")
        if "maximum" in schema and value > schema["maximum"]:
            raise ToolInputError(f"{path}: en fazla {schema['maximum']} olmalı")
    if isinstance(value, str):
        if "maxLength" in schema and len(value) > schema["maxLength"]:
            raise ToolInputError(f"{path}: en fazla {schema['maxLength']} karakter")
        if "minLength" in schema and len(value) < schema["minLength"]:
            raise ToolInputError(f"{path}: en az {schema['minLength']} karakter")
    if isinstance(value, list):
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            raise ToolInputError(f"{path}: en fazla {schema['maxItems']} öğe")
        if "minItems" in schema and len(value) < schema["minItems"]:
            raise ToolInputError(f"{path}: en az {schema['minItems']} öğe")
        if "items" in schema:
            for i, item in enumerate(value):
                validate(item, schema["items"], f"{path}[{i}]")
    if isinstance(value, dict):
        props = schema.get("properties", {})
        for key in schema.get("required", []):
            if key not in value:
                raise ToolInputError(f"{path}: '{key}' alanı gerekli")
        if schema.get("additionalProperties") is False:
            extra = set(value) - set(props)
            if extra:
                raise ToolInputError(f"{path}: tanınmayan alan(lar): {sorted(extra)}")
        for key, sub in props.items():
            if key in value:
                validate(value[key], sub, f"{path}.{key}")
    return value


def apply_defaults(value, schema):
    if isinstance(value, dict) and schema.get("type") == "object":
        for key, sub in schema.get("properties", {}).items():
            if key not in value and "default" in sub:
                value[key] = json.loads(json.dumps(sub["default"]))
            if key in value:
                apply_defaults(value[key], sub)
    return value


# ---------------------------------------------------------------------------
# Registry objects
# ---------------------------------------------------------------------------
class Tool:
    def __init__(self, name, title, description, input_schema, handler, *, read_only=True,
                 open_world=False, idempotent=True, output_schema=None):
        self.name, self.title, self.description = name, title, description
        self.input_schema, self.handler, self.output_schema = input_schema, handler, output_schema
        self.annotations = {"title": title, "readOnlyHint": read_only, "destructiveHint": False,
                            "idempotentHint": idempotent, "openWorldHint": open_world}

    def describe(self, version):
        d = {"name": self.name, "description": self.description, "inputSchema": self.input_schema}
        if version >= "2025-03-26":
            d["annotations"] = self.annotations
        if version >= "2025-06-18":
            d["title"] = self.title
            if self.output_schema:
                d["outputSchema"] = self.output_schema
        return d


class Resource:
    def __init__(self, uri, name, description, mime, reader):
        self.uri, self.name, self.description, self.mime, self.reader = uri, name, description, mime, reader


class Prompt:
    def __init__(self, name, description, arguments, render):
        self.name, self.description, self.arguments, self.render = name, description, arguments, render


class Registry:
    def __init__(self, name, title, version, instructions):
        self.name, self.title, self.version, self.instructions = name, title, version, instructions
        self.tools, self.resources, self.prompts = {}, {}, {}

    def tool(self, t):
        self.tools[t.name] = t

    def resource(self, r):
        self.resources[r.uri] = r

    def prompt(self, p):
        self.prompts[p.name] = p

    # -- tool execution shared by MCP and REST ------------------------------------------
    def call(self, name, arguments):
        tool = self.tools.get(name)
        if tool is None:
            raise KeyError(name)
        args = apply_defaults(dict(arguments or {}), tool.input_schema)
        validate(args, tool.input_schema)
        started = time.monotonic()
        result = tool.handler(args)
        result.setdefault("structured", {})
        result["structured"].setdefault("_meta", {})["elapsed_ms"] = round(1000 * (time.monotonic() - started))
        return result


# ---------------------------------------------------------------------------
# JSON-RPC session
# ---------------------------------------------------------------------------
class Session:
    def __init__(self, registry):
        self.registry = registry
        self.version = LATEST
        self.initialized = False

    def handle(self, msg):
        """Return a response dict, or None for notifications/responses."""
        if not isinstance(msg, dict) or msg.get("jsonrpc") != "2.0":
            return _error(None, INVALID_REQUEST, "Geçersiz JSON-RPC iletisi")
        if "method" not in msg:
            return None  # a response to a server->client request; we never send any
        mid, method, params = msg.get("id"), msg["method"], msg.get("params") or {}
        is_note = "id" not in msg
        try:
            result = self._dispatch(method, params)
        except _RpcError as exc:
            return None if is_note else _error(mid, exc.code, exc.message, exc.data)
        except Exception as exc:  # pragma: no cover - defensive
            traceback.print_exc(file=sys.stderr)
            return None if is_note else _error(mid, INTERNAL_ERROR, f"İç hata: {exc}")
        if is_note:
            return None
        return {"jsonrpc": "2.0", "id": mid, "result": result}

    def _dispatch(self, method, params):
        reg = self.registry
        if method == "initialize":
            asked = params.get("protocolVersion", LATEST)
            self.version = asked if asked in SUPPORTED_VERSIONS else LATEST
            self.initialized = True
            return {"protocolVersion": self.version,
                    "capabilities": {"tools": {"listChanged": False},
                                     "resources": {"subscribe": False, "listChanged": False},
                                     "prompts": {"listChanged": False}, "logging": {}},
                    "serverInfo": {"name": reg.name, "title": reg.title, "version": reg.version},
                    "instructions": reg.instructions}
        if method.startswith("notifications/"):
            return {}
        if method == "ping":
            return {}
        if method == "logging/setLevel":
            return {}
        if method == "tools/list":
            return {"tools": [t.describe(self.version) for t in reg.tools.values()]}
        if method == "tools/call":
            name = params.get("name")
            if name not in reg.tools:
                raise _RpcError(INVALID_PARAMS, f"Bilinmeyen araç: {name}")
            return self._call_tool(name, params.get("arguments") or {})
        if method == "resources/list":
            return {"resources": [{"uri": r.uri, "name": r.name, "description": r.description,
                                   "mimeType": r.mime} for r in reg.resources.values()]}
        if method == "resources/templates/list":
            return {"resourceTemplates": [
                {"uriTemplate": "oligoscout://rule/{rule_id}", "name": "QC kuralı",
                 "description": "Tek bir QC kuralının künyesi (ör. oligoscout://rule/TH-02)",
                 "mimeType": "application/json"}]}
        if method == "resources/read":
            uri = params.get("uri", "")
            r = reg.resources.get(uri)
            if r is None and uri.startswith("oligoscout://rule/"):
                r = reg.resources.get("oligoscout://catalog")
                text = r.reader(uri.rsplit("/", 1)[-1]) if r else None
                if text is None:
                    raise _RpcError(-32002, f"Kaynak bulunamadı: {uri}")
                return {"contents": [{"uri": uri, "mimeType": "application/json", "text": text}]}
            if r is None:
                raise _RpcError(-32002, f"Kaynak bulunamadı: {uri}")
            return {"contents": [{"uri": r.uri, "mimeType": r.mime, "text": r.reader(None)}]}
        if method == "prompts/list":
            return {"prompts": [{"name": p.name, "description": p.description, "arguments": p.arguments}
                                for p in reg.prompts.values()]}
        if method == "prompts/get":
            p = reg.prompts.get(params.get("name"))
            if p is None:
                raise _RpcError(INVALID_PARAMS, "Bilinmeyen istem")
            return p.render(params.get("arguments") or {})
        if method == "completion/complete":
            return {"completion": {"values": [], "hasMore": False}}
        raise _RpcError(METHOD_NOT_FOUND, f"Yöntem desteklenmiyor: {method}")

    def _call_tool(self, name, arguments):
        try:
            out = self.registry.call(name, arguments)
        except ToolInputError as exc:
            return {"content": [{"type": "text", "text": f"Girdi hatası: {exc}"}], "isError": True}
        except ValueError as exc:
            return {"content": [{"type": "text", "text": f"Analiz tamamlanamadı: {exc}"}], "isError": True}
        except Exception as exc:  # pragma: no cover
            traceback.print_exc(file=sys.stderr)
            return {"content": [{"type": "text", "text": f"Motor hatası: {exc}"}], "isError": True}
        structured = out.get("structured", {})
        text = out.get("text") or json.dumps(structured, ensure_ascii=False)[:4000]
        result = {"content": [{"type": "text", "text": text}], "isError": bool(out.get("error"))}
        if self.version >= "2025-06-18":
            result["structuredContent"] = structured
        else:  # older clients: also ship the JSON as a second text block
            result["content"].append({"type": "text", "text": json.dumps(structured, ensure_ascii=False)})
        return result


class _RpcError(Exception):
    def __init__(self, code, message, data=None):
        super().__init__(message)
        self.code, self.message, self.data = code, message, data


def _error(mid, code, message, data=None):
    err = {"code": code, "message": message}
    if data is not None:
        err["data"] = data
    return {"jsonrpc": "2.0", "id": mid, "error": err}


# ---------------------------------------------------------------------------
# stdio transport
# ---------------------------------------------------------------------------
def serve_stdio(registry):
    session = Session(registry)
    out = sys.stdout
    lock = threading.Lock()
    log(f"stdio hazır · {len(registry.tools)} araç")
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            resp = _error(None, PARSE_ERROR, "JSON ayrıştırılamadı")
        else:
            if isinstance(msg, list):  # 2025-03-26 batch support
                resp = [r for r in (session.handle(m) for m in msg) if r is not None] or None
            else:
                resp = session.handle(msg)
        if resp is not None:
            with lock:
                out.write(json.dumps(resp, ensure_ascii=False) + "\n")
                out.flush()


# ---------------------------------------------------------------------------
# Streamable HTTP transport + REST/OpenAPI
# ---------------------------------------------------------------------------
class HttpServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, addr, registry, *, token="", allowed_origins=(), public_url="", rest=True, oauth=None):
        super().__init__(addr, _Handler)
        self.registry, self.token, self.allowed_origins = registry, token, allowed_origins
        self.oauth = oauth
        self.public_url, self.rest = public_url.rstrip("/"), rest
        self.sessions, self.lock = {}, threading.Lock()
        self.hits = {}


class _Handler(BaseHTTPRequestHandler):
    server_version = "OligoScoutMCP"
    protocol_version = "HTTP/1.1"
    RATE = (60, 60.0)  # 60 calls per 60 s per client address

    def log_message(self, fmt, *args):
        log(self.address_string(), fmt % args)

    def _send(self, status, payload=None, ctype="application/json", headers=None):
        body = b"" if payload is None else (
            payload if isinstance(payload, bytes) else json.dumps(payload, ensure_ascii=False).encode())
        self.send_response(status)
        if body:
            self.send_header("Content-Type", ctype + ("; charset=utf-8" if "json" in ctype or "text" in ctype else ""))
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        for k, v in (headers or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if body:
            self.wfile.write(body)

    def _guard(self):
        s = self.server
        origin = (self.headers.get("Origin") or "").rstrip("/")
        if origin and origin != "null":
            host = urlsplit(origin).hostname or ""
            local = host in ("localhost", "127.0.0.1", "::1")
            if origin not in s.allowed_origins and not local:
                # DNS-rebinding guard required by the MCP Streamable HTTP spec
                self._send(HTTPStatus.FORBIDDEN, {"error": "Origin izinli değil"})
                return False
        if s.token or s.oauth:
            auth = self.headers.get("Authorization", "")
            given = auth[7:].strip() if auth.lower().startswith("bearer ") else self.headers.get("X-API-Key", "")
            ok = bool(s.token) and hmac.compare_digest(given.encode(), s.token.encode())
            ok = ok or (s.oauth is not None and s.oauth.valid(given))
            if not ok:
                challenge = s.oauth.www_authenticate() if s.oauth else 'Bearer realm="oligoscout"'
                self._send(HTTPStatus.UNAUTHORIZED, {"error": "Kimlik doğrulama gerekli"},
                           headers={"WWW-Authenticate": challenge})
                return False
        return self._rate_ok()

    def _rate_ok(self):
        s = self.server
        ip, now = self.client_address[0], time.monotonic()
        with s.lock:
            window = [t for t in s.hits.get(ip, []) if now - t < self.RATE[1]]
            window.append(now)
            s.hits[ip] = window
        if len(window) > self.RATE[0]:
            self._send(HTTPStatus.TOO_MANY_REQUESTS, {"error": "Çok fazla istek"}, headers={"Retry-After": "30"})
            return False
        return True

    def _body(self, limit=256_000):
        size = int(self.headers.get("Content-Length") or 0)
        if not 0 < size <= limit:
            raise ValueError("Gövde boş ya da çok büyük")
        return json.loads(self.rfile.read(size))

    # -- OAuth 2.1 (mcp_server/oauth.py) --------------------------------------
    _CORS = {"Access-Control-Allow-Origin": "*"}

    def _oauth_get(self, path, query):
        o = self.server.oauth
        if path in ("/.well-known/oauth-protected-resource", "/.well-known/oauth-protected-resource/mcp",
                    "/mcp/.well-known/oauth-protected-resource"):
            self._send(HTTPStatus.OK, o.protected_resource(), headers=self._CORS)
            return True
        if path in ("/.well-known/oauth-authorization-server", "/.well-known/oauth-authorization-server/mcp",
                    "/.well-known/openid-configuration", "/mcp/.well-known/oauth-authorization-server"):
            self._send(HTTPStatus.OK, o.authorization_server(), headers=self._CORS)
            return True
        if path.endswith("/oauth/authorize"):
            q = {k: v[0] for k, v in parse_qs(query).items()}
            err, client = o.authorize_check(q)
            if err:
                self._send(HTTPStatus.BAD_REQUEST, f"OAuth isteği geçersiz: {err}".encode(), ctype="text/plain")
            else:
                # form-action also governs the redirect that follows the form POST: browsers
                # block the 302 to the client's callback unless its origin is listed here.
                cb = urlsplit(q.get("redirect_uri", ""))
                cb_src = (f"{cb.scheme}://{cb.netloc}" if cb.scheme in ("http", "https") and cb.netloc
                          else f"{cb.scheme}:" if cb.scheme else "")
                self._send(HTTPStatus.OK, o.consent_page(q, client).encode(), ctype="text/html",
                           headers={"Content-Security-Policy": "default-src 'none'; style-src 'unsafe-inline'; "
                                    f"form-action 'self' {cb_src}; frame-ancestors 'none'"})
            return True
        return False

    def _oauth_post(self, path):
        o = self.server.oauth
        if not path.endswith(("/oauth/register", "/oauth/authorize", "/oauth/token")):
            return False
        if not self._rate_ok():
            return True
        size = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(size) if 0 < size <= 64_000 else b""
        try:
            form = parse_form(raw, self.headers.get("Content-Type", ""))
        except (ValueError, json.JSONDecodeError):
            self._send(HTTPStatus.BAD_REQUEST, {"error": "invalid_request"}, headers=self._CORS)
            return True
        if path.endswith("/oauth/register"):
            status, body = o.register(json.loads(raw or b"{}") if raw else {})
            self._send(status, body, headers=self._CORS)
        elif path.endswith("/oauth/token"):
            status, body = o.token(form)
            self._send(status, body, headers={**self._CORS, "Pragma": "no-cache"})
        else:
            status, where = o.authorize_decide(form)
            if status == 302:
                self._send(HTTPStatus.FOUND, None, headers={"Location": where})
            else:
                self._send(status, f"OAuth: {where}".encode(), ctype="text/plain")
        return True

    def do_GET(self):
        parts = urlsplit(self.path)
        path = parts.path
        if self.server.oauth and self._oauth_get(path, parts.query):
            return
        if path in ("/healthz", "/mcp/healthz"):
            reg = self.server.registry
            return self._send(HTTPStatus.OK, {"status": "ok", "server": reg.name, "version": reg.version,
                                              "tools": len(reg.tools)})
        if path in ("/openapi.json", "/mcp/openapi.json") and self.server.rest:
            return self._send(HTTPStatus.OK, openapi(self.server.registry, self.server.public_url))
        if path in ("/mcp", "/mcp/"):
            # We do not offer a server-initiated SSE stream (spec permits 405).
            return self._send(HTTPStatus.METHOD_NOT_ALLOWED, None, headers={"Allow": "POST, DELETE"})
        self._send(HTTPStatus.NOT_FOUND, {"error": "Bulunamadı"})

    def do_DELETE(self):
        if not self._guard():
            return
        sid = self.headers.get("Mcp-Session-Id")
        with self.server.lock:
            self.server.sessions.pop(sid, None)
        self._send(HTTPStatus.OK, None)

    def do_OPTIONS(self):
        self._send(HTTPStatus.NO_CONTENT, None, headers={
            "Access-Control-Allow-Methods": "POST, GET, DELETE, OPTIONS",
            "Access-Control-Allow-Headers": "Content-Type, Authorization, Mcp-Session-Id, MCP-Protocol-Version"})

    def do_POST(self):
        path = urlsplit(self.path).path
        if self.server.oauth and self._oauth_post(path):
            return
        if not self._guard():
            return
        try:
            body = self._body()
        except (ValueError, json.JSONDecodeError) as exc:
            return self._send(HTTPStatus.BAD_REQUEST, _error(None, PARSE_ERROR, str(exc)))
        if path.startswith(("/api/tools/", "/mcp/api/tools/")) and self.server.rest:
            return self._rest(path.rsplit("/", 1)[-1], body)
        if path not in ("/mcp", "/mcp/"):
            return self._send(HTTPStatus.NOT_FOUND, {"error": "Bulunamadı"})
        s = self.server
        sid = self.headers.get("Mcp-Session-Id")
        with s.lock:
            session = s.sessions.get(sid) if sid else None
        msgs = body if isinstance(body, list) else [body]
        headers = {}
        if session is None:
            session = Session(s.registry)
            if any(isinstance(m, dict) and m.get("method") == "initialize" for m in msgs):
                sid = secrets.token_hex(16)
                with s.lock:
                    s.sessions[sid] = session
                    if len(s.sessions) > 2000:
                        s.sessions.pop(next(iter(s.sessions)))
                headers["Mcp-Session-Id"] = sid
            # Stateless fallback: requests without a known session are still served.
        responses = [r for r in (session.handle(m) for m in msgs) if r is not None]
        if not responses:
            return self._send(HTTPStatus.ACCEPTED, None, headers=headers)
        payload = responses if isinstance(body, list) else responses[0]
        self._send(HTTPStatus.OK, payload, headers=headers)

    def _rest(self, name, body):
        reg = self.server.registry
        if name not in reg.tools:
            return self._send(HTTPStatus.NOT_FOUND, {"error": f"Bilinmeyen araç: {name}"})
        try:
            out = reg.call(name, body if isinstance(body, dict) else {})
        except (ToolInputError, ValueError) as exc:
            return self._send(HTTPStatus.UNPROCESSABLE_ENTITY, {"error": str(exc)})
        except Exception as exc:  # pragma: no cover
            traceback.print_exc(file=sys.stderr)
            return self._send(HTTPStatus.INTERNAL_SERVER_ERROR, {"error": f"Motor hatası: {exc}"})
        self._send(HTTPStatus.OK, {"summary": out.get("text", ""), "result": out.get("structured", {})})


def serve_http(registry, host="127.0.0.1", port=8766, **kw):
    httpd = HttpServer((host, port), registry, **kw)
    log(f"HTTP hazır · http://{host}:{port}/mcp · REST /api/tools/<ad> · /openapi.json")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


def openapi(registry, public_url=""):
    paths = {}
    for t in registry.tools.values():
        paths[f"/api/tools/{t.name}"] = {"post": {
            "operationId": t.name, "summary": t.title,
            "description": t.description[:300],
            "x-openai-isConsequential": False,
            "requestBody": {"required": True, "content": {"application/json": {"schema": t.input_schema}}},
            "responses": {"200": {"description": "Kanıt kayıtlı sonuç",
                                  "content": {"application/json": {"schema": {
                                      "type": "object", "properties": {
                                          "summary": {"type": "string"},
                                          "result": {"type": "object", "additionalProperties": True}}}}}},
                          "422": {"description": "Girdi hatası"}}}}
    return {"openapi": "3.1.0",
            "info": {"title": registry.title, "version": registry.version,
                     "description": registry.instructions[:1000]},
            "servers": [{"url": public_url or os.getenv("OLIGOSCOUT_PUBLIC_URL", "http://127.0.0.1:8766")}],
            "paths": paths,
            "components": {"securitySchemes": {"bearer": {"type": "http", "scheme": "bearer"}}},
            "security": [{"bearer": []}]}
