#!/bin/sh
# OligoScout container entrypoint: `mcp` (default), `web`, or any command.
set -e
cd /opt/oligoscout
case "${1:-mcp}" in
  mcp)
    shift 2>/dev/null || true
    exec python3 mcp_server/server.py --http --host 0.0.0.0 --port "${OLIGOSCOUT_MCP_PORT:-8766}" "$@"
    ;;
  web)
    shift
    exec python3 web/local_server.py --host 0.0.0.0 --port "${OLIGOSCOUT_WEB_PORT:-8765}" "$@"
    ;;
  *)
    exec "$@"
    ;;
esac
