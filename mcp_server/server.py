#!/usr/bin/env python3
"""OligoScout MCP server.

    python mcp_server/server.py                     # stdio (Claude Desktop/Code, Cursor, VS Code, Gemini CLI ...)
    python mcp_server/server.py --http --port 8766  # Streamable HTTP + REST + OpenAPI (ChatGPT, claude.ai, remote)

Environment:
    OLIGOSCOUT_ROOT             project root (default: parent of this folder)
    OLIGOSCOUT_PANEL_BINARY     path to oligoscout_panel if not in OligoScout/build/release
    OLIGOSCOUT_MCP_TOKEN        require "Authorization: Bearer <token>" on HTTP
    OLIGOSCOUT_ALLOWED_ORIGINS  comma-separated browser origins allowed on HTTP (DNS-rebinding guard)
    OLIGOSCOUT_PUBLIC_URL       public base URL advertised in /openapi.json
    OLIGOSCOUT_GNOMAD_API=0     offline mode: only the local gnomAD cache is used
    OLIGOSCOUT_OAUTH=1          enable OAuth 2.1 (PKCE + dynamic client registration) for public clients
    OLIGOSCOUT_OAUTH_INVITE     optional invite code required on the consent page
    OLIGOSCOUT_OAUTH_STATE      JSON file that keeps clients/tokens across restarts (tokens stored hashed)
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from protocol import serve_http, serve_stdio  # noqa: E402
from tools import build_registry  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="OligoScout MCP server")
    ap.add_argument("--http", action="store_true", help="serve Streamable HTTP instead of stdio")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8766)
    ap.add_argument("--no-rest", action="store_true", help="disable /api/tools and /openapi.json")
    args = ap.parse_args()
    reg = build_registry()
    if args.http:
        origins = tuple(o.strip().rstrip("/") for o in os.getenv("OLIGOSCOUT_ALLOWED_ORIGINS", "").split(",") if o.strip())
        if args.host not in ("127.0.0.1", "localhost", "::1") and not (os.getenv("OLIGOSCOUT_MCP_TOKEN")
                                                                          or os.getenv("OLIGOSCOUT_OAUTH") == "1"):
            print("[oligoscout-mcp] UYARI: herkese açık arayüzde token yok; OLIGOSCOUT_MCP_TOKEN ayarlayın "
                  "ya da yalnızca TLS sonlandıran bir ters vekil arkasında çalıştırın.", file=sys.stderr)
        oauth = None
        if os.getenv("OLIGOSCOUT_OAUTH") == "1":
            from oauth import OAuthServer
            oauth = OAuthServer(os.getenv("OLIGOSCOUT_PUBLIC_URL", ""), invite=os.getenv("OLIGOSCOUT_OAUTH_INVITE", ""),
                                state_path=os.getenv("OLIGOSCOUT_OAUTH_STATE", ""))
            print(f"[oligoscout-mcp] OAuth 2.1 açık · issuer {oauth.issuer} · kaynak {oauth.resource}", file=sys.stderr)
        serve_http(reg, args.host, args.port, token=os.getenv("OLIGOSCOUT_MCP_TOKEN", ""),
                   allowed_origins=origins, public_url=os.getenv("OLIGOSCOUT_PUBLIC_URL", ""), rest=not args.no_rest,
                   oauth=oauth)
    else:
        serve_stdio(reg)


if __name__ == "__main__":
    main()
