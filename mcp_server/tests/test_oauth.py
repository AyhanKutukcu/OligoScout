"""End-to-end OAuth 2.1 flow against a live HTTP server (no network beyond loopback).

    python3 -m unittest discover -s mcp_server/tests
"""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import secrets
import sys
import threading
import unittest
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from oauth import OAuthServer  # noqa: E402
from protocol import HttpServer, Registry, Tool  # noqa: E402


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k):
        return None


OPENER = urllib.request.build_opener(_NoRedirect)


def req(url, data=None, headers=None, method=None):
    r = urllib.request.Request(url, data=data, headers=headers or {}, method=method)
    try:
        with OPENER.open(r, timeout=10) as resp:
            return resp.status, dict(resp.headers), resp.read()
    except urllib.error.HTTPError as e:
        return e.code, dict(e.headers), e.read()


class OAuthFlow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        reg = Registry("oligoscout-test", "0", "OligoScout test", "test")
        reg.tool(Tool("echo", "Echo", "echo", {"type": "object", "properties": {}}, lambda a: {"text": "ok"}))
        probe = HttpServer(("127.0.0.1", 0), reg)
        port = probe.server_address[1]
        probe.server_close()
        cls.base = f"http://127.0.0.1:{port}"
        oauth = OAuthServer(cls.base + "/mcp")
        cls.httpd = HttpServer(("127.0.0.1", port), reg, oauth=oauth, public_url=cls.base + "/mcp")
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()

    def init_body(self):
        return json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
            "protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "t", "version": "1"}}}).encode()

    def test_full_flow(self):
        st, h, _ = req(self.base + "/mcp", self.init_body(), {"Content-Type": "application/json"})
        self.assertEqual(st, 401)
        meta_url = re.search(r'resource_metadata="([^"]+)"', h["WWW-Authenticate"]).group(1)
        st, _, body = req(meta_url)
        self.assertEqual(st, 200)
        prm = json.loads(body)
        self.assertEqual(prm["resource"], self.base + "/mcp")
        st, _, body = req(prm["authorization_servers"][0] + "/.well-known/oauth-authorization-server")
        asm = json.loads(body)
        self.assertIn("S256", asm["code_challenge_methods_supported"])
        redirect = "http://127.0.0.1:9/callback"
        st, _, body = req(asm["registration_endpoint"], json.dumps({"redirect_uris": [redirect],
                                                                    "client_name": "Test"}).encode(),
                          {"Content-Type": "application/json"})
        self.assertEqual(st, 201)
        cid = json.loads(body)["client_id"]
        verifier = secrets.token_urlsafe(48)
        challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).rstrip(b"=").decode()
        q = {"response_type": "code", "client_id": cid, "redirect_uri": redirect, "code_challenge": challenge,
             "code_challenge_method": "S256", "state": "xyz", "scope": "oligoscout", "resource": self.base + "/mcp"}
        st, _, page = req(asm["authorization_endpoint"] + "?" + urllib.parse.urlencode(q))
        self.assertEqual(st, 200)
        fields = dict(re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', page.decode()))
        fields = {k: v.replace("&amp;", "&") for k, v in fields.items()}
        fields["decision"] = "allow"
        st, h, _ = req(asm["authorization_endpoint"], urllib.parse.urlencode(fields).encode(),
                       {"Content-Type": "application/x-www-form-urlencoded"})
        self.assertEqual(st, 302)
        loc = urllib.parse.urlsplit(h["Location"])
        params = dict(urllib.parse.parse_qsl(loc.query))
        self.assertEqual(params["state"], "xyz")
        bad = {"grant_type": "authorization_code", "code": params["code"], "redirect_uri": redirect,
               "client_id": cid, "code_verifier": "x" * 50}
        st, _, body = req(asm["token_endpoint"], urllib.parse.urlencode(bad).encode(),
                          {"Content-Type": "application/x-www-form-urlencoded"})
        self.assertEqual(st, 400)                       # wrong verifier burns the code
        # a fresh code with the right verifier
        st, _, page = req(asm["authorization_endpoint"] + "?" + urllib.parse.urlencode(q))
        fields = dict(re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', page.decode()))
        fields = {k: v.replace("&amp;", "&") for k, v in fields.items()}
        fields["decision"] = "allow"
        _, h, _ = req(asm["authorization_endpoint"], urllib.parse.urlencode(fields).encode(),
                      {"Content-Type": "application/x-www-form-urlencoded"})
        code = dict(urllib.parse.parse_qsl(urllib.parse.urlsplit(h["Location"]).query))["code"]
        good = {**bad, "code": code, "code_verifier": verifier}
        st, _, body = req(asm["token_endpoint"], urllib.parse.urlencode(good).encode(),
                          {"Content-Type": "application/x-www-form-urlencoded"})
        self.assertEqual(st, 200)
        tok = json.loads(body)
        st, _, body = req(self.base + "/mcp", self.init_body(),
                          {"Content-Type": "application/json", "Authorization": "Bearer " + tok["access_token"]})
        self.assertEqual(st, 200)
        self.assertEqual(json.loads(body)["result"]["serverInfo"]["name"], "oligoscout-test")
        # refresh rotation: the old refresh token works once
        ref = {"grant_type": "refresh_token", "refresh_token": tok["refresh_token"], "client_id": cid}
        st, _, body = req(asm["token_endpoint"], urllib.parse.urlencode(ref).encode(),
                          {"Content-Type": "application/x-www-form-urlencoded"})
        self.assertEqual(st, 200)
        st, _, _ = req(asm["token_endpoint"], urllib.parse.urlencode(ref).encode(),
                       {"Content-Type": "application/x-www-form-urlencoded"})
        self.assertEqual(st, 400)

    def test_unregistered_redirect_rejected(self):
        st, _, _ = req(self.base + "/mcp/oauth/register", json.dumps({"redirect_uris": ["http://evil.example/cb"]}).encode(),
                       {"Content-Type": "application/json"})
        self.assertEqual(st, 400)

    def test_deny_redirects_with_error(self):
        redirect = "http://localhost:9/cb"
        _, _, body = req(self.base + "/mcp/oauth/register", json.dumps({"redirect_uris": [redirect]}).encode(),
                         {"Content-Type": "application/json"})
        cid = json.loads(body)["client_id"]
        q = {"response_type": "code", "client_id": cid, "redirect_uri": redirect, "code_challenge": "a" * 43,
             "code_challenge_method": "S256", "state": "s1"}
        _, hdr, page = req(self.base + "/mcp/oauth/authorize?" + urllib.parse.urlencode(q))
        # the consent form's POST is answered by a redirect to the client: CSP form-action
        # must allow that origin or browsers silently block the "allow" button
        self.assertIn("form-action 'self' http://localhost:9;", hdr["Content-Security-Policy"])
        fields = dict(re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', page.decode()))
        fields = {k: v.replace("&amp;", "&") for k, v in fields.items()}
        fields["decision"] = "deny"
        st, h, _ = req(self.base + "/mcp/oauth/authorize", urllib.parse.urlencode(fields).encode(),
                       {"Content-Type": "application/x-www-form-urlencoded"})
        self.assertEqual(st, 302)
        self.assertIn("error=access_denied", h["Location"])


if __name__ == "__main__":
    unittest.main()
