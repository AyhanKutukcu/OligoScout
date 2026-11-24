"""Minimal OAuth 2.1 authorization server for the public OligoScout MCP endpoint.

Implements what remote MCP clients (ChatGPT connectors, claude.ai custom
connectors, MCP Inspector) need, with the standard library only:

  * RFC 9728  protected-resource metadata   /.well-known/oauth-protected-resource[/mcp]
  * RFC 8414  authorization-server metadata /.well-known/oauth-authorization-server
  * RFC 7591  dynamic client registration   POST <prefix>/oauth/register   (public clients)
  * OAuth 2.1 authorization code + PKCE S256 GET/POST <prefix>/oauth/authorize (consent page)
              token endpoint                 POST <prefix>/oauth/token  (authorization_code, refresh_token)
  * Bearer validation for /mcp and the REST tools, with `WWW-Authenticate:
    Bearer resource_metadata="..."` on 401 so clients can discover the server.

Policy ("open registration", the default): anyone may register a client and
approve access on the consent page; every tool is read-only and the MCP server
rate-limits each address. Set OLIGOSCOUT_OAUTH_INVITE to require an invite code
on the consent page instead.

Tokens are opaque random strings; only their SHA-256 is kept. State survives a
restart when OLIGOSCOUT_OAUTH_STATE names a writable JSON file.
"""
from __future__ import annotations

import base64
import hashlib
import hmac
import html
import json
import os
import secrets
import threading
import time
from urllib.parse import parse_qs, urlencode, urlsplit

ACCESS_TTL = 3600
REFRESH_TTL = 30 * 24 * 3600
CODE_TTL = 300
SCOPE = "oligoscout"
MAX_CLIENTS = 5000


def _h(token: str) -> str:
    return hashlib.sha256(token.encode()).hexdigest()


def _b64url(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).rstrip(b"=").decode()


def _redirect_ok(uri: str) -> bool:
    try:
        u = urlsplit(uri)
    except ValueError:
        return False
    if u.fragment:
        return False
    if u.scheme == "https" and u.hostname:
        return True
    # loopback redirects for desktop/CLI clients (RFC 8252)
    return u.scheme == "http" and u.hostname in ("localhost", "127.0.0.1", "::1")


class OAuthServer:
    def __init__(self, public_url: str, *, issuer: str = "", invite: str = "", state_path: str = ""):
        pub = urlsplit(public_url)
        if not (pub.scheme and pub.netloc):
            raise ValueError("OAuth requires OLIGOSCOUT_PUBLIC_URL, e.g. https://oligoscout.org/mcp")
        self.resource = public_url.rstrip("/")                         # https://host/mcp
        self.issuer = (issuer or f"{pub.scheme}://{pub.netloc}").rstrip("/")
        prefix = pub.path.rstrip("/")                                  # /mcp
        self.base = f"{pub.scheme}://{pub.netloc}{prefix}"
        self.prefix = prefix
        self.invite = invite
        self.state_path = state_path
        self.lock = threading.Lock()
        self.secret = secrets.token_bytes(32)
        self.clients: dict[str, dict] = {}
        self.codes: dict[str, dict] = {}
        self.tokens: dict[str, dict] = {}      # sha256(token) -> record
        self._load()

    # -- persistence ---------------------------------------------------------
    def _load(self):
        if self.state_path and os.path.isfile(self.state_path):
            with open(self.state_path, encoding="utf-8") as fh:
                d = json.load(fh)
            self.clients, self.tokens = d.get("clients", {}), d.get("tokens", {})
            self.secret = bytes.fromhex(d.get("secret", self.secret.hex()))

    def _save(self):
        if not self.state_path:
            return
        now = time.time()
        self.tokens = {k: v for k, v in self.tokens.items() if v["exp"] > now}
        tmp = self.state_path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump({"clients": self.clients, "tokens": self.tokens, "secret": self.secret.hex()}, fh)
        os.replace(tmp, self.state_path)

    # -- metadata ------------------------------------------------------------
    def resource_metadata_url(self) -> str:
        return f"{self.issuer}/.well-known/oauth-protected-resource{self.prefix}"

    def protected_resource(self) -> dict:
        return {"resource": self.resource, "authorization_servers": [self.issuer],
                "scopes_supported": [SCOPE], "bearer_methods_supported": ["header"],
                "resource_name": "OligoScout primer design & QC (MCP)"}

    def authorization_server(self) -> dict:
        return {"issuer": self.issuer,
                "authorization_endpoint": f"{self.base}/oauth/authorize",
                "token_endpoint": f"{self.base}/oauth/token",
                "registration_endpoint": f"{self.base}/oauth/register",
                "response_types_supported": ["code"],
                "grant_types_supported": ["authorization_code", "refresh_token"],
                "code_challenge_methods_supported": ["S256"],
                "token_endpoint_auth_methods_supported": ["none"],
                "scopes_supported": [SCOPE]}

    def www_authenticate(self) -> str:
        return f'Bearer resource_metadata="{self.resource_metadata_url()}", scope="{SCOPE}"'

    # -- registration (RFC 7591) ----------------------------------------------
    def register(self, body: dict) -> tuple[int, dict]:
        uris = body.get("redirect_uris")
        if not isinstance(uris, list) or not uris or not all(isinstance(u, str) and _redirect_ok(u) for u in uris):
            return 400, {"error": "invalid_redirect_uri",
                         "error_description": "redirect_uris must be https or loopback http URLs"}
        name = str(body.get("client_name") or "MCP client")[:80]
        with self.lock:
            if len(self.clients) >= MAX_CLIENTS:
                oldest = min(self.clients, key=lambda k: self.clients[k]["created"])
                self.clients.pop(oldest)
            cid = "osc_" + secrets.token_urlsafe(18)
            self.clients[cid] = {"name": name, "redirect_uris": uris, "created": int(time.time())}
            self._save()
        return 201, {"client_id": cid, "client_name": name, "redirect_uris": uris,
                     "grant_types": ["authorization_code", "refresh_token"], "response_types": ["code"],
                     "token_endpoint_auth_method": "none", "client_id_issued_at": int(time.time())}

    # -- authorization ---------------------------------------------------------
    def _csrf(self, params: str) -> str:
        return hmac.new(self.secret, params.encode(), hashlib.sha256).hexdigest()

    def authorize_check(self, q: dict) -> tuple[str | None, dict]:
        """Validate an authorization request. Returns (error, client)."""
        cid = q.get("client_id", "")
        client = self.clients.get(cid)
        if not client:
            return "unknown client_id", {}
        if q.get("redirect_uri") not in client["redirect_uris"]:
            return "redirect_uri not registered for this client", {}
        if q.get("response_type") != "code":
            return "response_type must be code", client
        if q.get("code_challenge_method") != "S256" or len(q.get("code_challenge", "")) < 43:
            return "PKCE S256 code_challenge required", client
        res = q.get("resource")
        if res and res.rstrip("/") != self.resource:
            return f"resource must be {self.resource}", client
        return None, client

    def consent_page(self, q: dict, client: dict) -> str:
        keep = {k: q.get(k, "") for k in ("client_id", "redirect_uri", "response_type", "code_challenge",
                                          "code_challenge_method", "state", "scope", "resource")}
        packed = urlencode(keep)
        hidden = "".join(f'<input type="hidden" name="{k}" value="{html.escape(v)}">' for k, v in keep.items())
        invite = ('<label>Davet kodu / invite code<input name="invite" autocomplete="off" required></label>'
                  if self.invite else "")
        host = html.escape(urlsplit(q.get("redirect_uri", "")).netloc)
        return f"""<!doctype html><html lang="tr"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>OligoScout · erişim izni</title>
<style>body{{font-family:system-ui,sans-serif;background:#f4f6f5;color:#1c2421;margin:0;padding:24px}}
main{{max-width:460px;margin:8vh auto;background:#fff;border:1px solid #d9e0dc;border-radius:12px;padding:28px}}
h1{{font-size:20px;margin:0 0 8px}}p{{line-height:1.55}}code{{background:#eef2f0;padding:1px 5px;border-radius:4px}}
label{{display:grid;gap:6px;margin:14px 0}}input{{padding:9px;border:1px solid #c4ccc8;border-radius:8px}}
.row{{display:flex;gap:10px;margin-top:18px}}button{{flex:1;padding:11px;border-radius:8px;border:1px solid #1f6f5c;font-size:15px;cursor:pointer}}
.ok{{background:#1f6f5c;color:#fff}}.no{{background:#fff;color:#1f6f5c}}small{{color:#5b6863}}</style>
<main><h1>OligoScout'a bağlan</h1>
<p><b>{html.escape(client["name"])}</b> (<code>{host}</code>) OligoScout primer tasarım ve QC araçlarını
senin adına kullanmak istiyor. Araçlar yalnızca okur; hesabına ya da dosyalarına erişmez.</p>
<p><small>{html.escape(client["name"])} wants to use OligoScout's read-only primer design and QC tools.
Results are computational (REP-03: not wet-lab calibrated).</small></p>
<form method="post" action="{html.escape(self.base)}/oauth/authorize">{hidden}
<input type="hidden" name="csrf" value="{self._csrf(packed)}">{invite}
<div class="row"><button class="no" name="decision" value="deny">Reddet</button>
<button class="ok" name="decision" value="allow">İzin ver</button></div></form></main></html>"""

    def authorize_decide(self, form: dict) -> tuple[int, str]:
        """Returns (status, redirect location or error text)."""
        keep = {k: form.get(k, "") for k in ("client_id", "redirect_uri", "response_type", "code_challenge",
                                             "code_challenge_method", "state", "scope", "resource")}
        if not hmac.compare_digest(form.get("csrf", ""), self._csrf(urlencode(keep))):
            return 400, "invalid form"
        err, _ = self.authorize_check(keep)
        if err:
            return 400, err
        sep = "&" if "?" in keep["redirect_uri"] else "?"
        if form.get("decision") != "allow":
            return 302, keep["redirect_uri"] + sep + urlencode({"error": "access_denied", "state": keep["state"]})
        if self.invite and not hmac.compare_digest(form.get("invite", ""), self.invite):
            return 403, "invalid invite code"
        code = secrets.token_urlsafe(32)
        with self.lock:
            self.codes[_h(code)] = {"client_id": keep["client_id"], "redirect_uri": keep["redirect_uri"],
                                    "challenge": keep["code_challenge"], "exp": time.time() + CODE_TTL,
                                    "scope": SCOPE, "resource": self.resource}
        return 302, keep["redirect_uri"] + sep + urlencode({"code": code, "state": keep["state"]})

    # -- token endpoint ---------------------------------------------------------
    def _issue(self, client_id: str) -> dict:
        access, refresh = secrets.token_urlsafe(32), secrets.token_urlsafe(32)
        now = time.time()
        self.tokens[_h(access)] = {"kind": "access", "client_id": client_id, "exp": now + ACCESS_TTL}
        self.tokens[_h(refresh)] = {"kind": "refresh", "client_id": client_id, "exp": now + REFRESH_TTL}
        self._save()
        return {"access_token": access, "token_type": "Bearer", "expires_in": ACCESS_TTL,
                "refresh_token": refresh, "scope": SCOPE}

    def token(self, form: dict) -> tuple[int, dict]:
        grant = form.get("grant_type")
        with self.lock:
            if grant == "authorization_code":
                rec = self.codes.pop(_h(form.get("code", "")), None)
                if not rec or rec["exp"] < time.time():
                    return 400, {"error": "invalid_grant"}
                if form.get("client_id") != rec["client_id"] or form.get("redirect_uri") != rec["redirect_uri"]:
                    return 400, {"error": "invalid_grant"}
                verifier = form.get("code_verifier", "")
                if not 43 <= len(verifier) <= 128 or \
                        _b64url(hashlib.sha256(verifier.encode()).digest()) != rec["challenge"]:
                    return 400, {"error": "invalid_grant", "error_description": "PKCE verification failed"}
                return 200, self._issue(rec["client_id"])
            if grant == "refresh_token":
                key = _h(form.get("refresh_token", ""))
                rec = self.tokens.get(key)
                if not rec or rec["kind"] != "refresh" or rec["exp"] < time.time() \
                        or form.get("client_id", rec["client_id"]) != rec["client_id"]:
                    return 400, {"error": "invalid_grant"}
                self.tokens.pop(key)                 # rotation: a refresh token is single-use
                return 200, self._issue(rec["client_id"])
        return 400, {"error": "unsupported_grant_type"}

    def valid(self, bearer: str) -> bool:
        rec = self.tokens.get(_h(bearer)) if bearer else None
        return bool(rec and rec["kind"] == "access" and rec["exp"] > time.time())


def parse_form(raw: bytes, ctype: str) -> dict:
    if "json" in (ctype or ""):
        d = json.loads(raw or b"{}")
        return {k: str(v) for k, v in d.items()} if isinstance(d, dict) else {}
    return {k: v[0] for k, v in parse_qs(raw.decode("utf-8", "replace")).items()}
