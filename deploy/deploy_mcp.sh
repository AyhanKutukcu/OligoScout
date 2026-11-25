#!/usr/bin/env bash
# Deploy / update the OligoScout MCP server on the oligoscout.org VPS.
# Run ON THE SERVER as a sudo-capable user, from the repository checkout:
#     sudo bash deploy/deploy_mcp.sh            # default checkout /opt/oligoscout
# It does not handle credentials: an optional static token lives in
# /etc/oligoscout/mcp.env, OAuth state in /var/lib/oligoscout/.
set -euo pipefail
APP=${APP:-/opt/oligoscout}
cd "$APP"
echo "== build C++ tools"
( cd OligoScout && cmake --preset release >/dev/null && \
  cmake --build --preset release --target oligoscout_panel oligoscout_contig_scan oligoscout_vpcr_scan -j"$(nproc)" )
echo "== python tests"
( cd reference_py && python3 -B -m unittest discover -s tests -q )
python3 -B -m unittest discover -s mcp_server/tests -q
echo "== extra contigs for SPC-08 (once)"
[[ -s OligoScout/data/annotations/grch38_extra_contigs.fa ]] || python3 reference_py/cli/build_extra_contigs.py
echo "== systemd"
install -d -m 750 -o oligoscout -g oligoscout /etc/oligoscout
[[ -f /etc/oligoscout/mcp.env ]] || install -m 600 -o oligoscout -g oligoscout /dev/null /etc/oligoscout/mcp.env
cp deploy/oligoscout-mcp.service /etc/systemd/system/oligoscout-mcp.service
systemctl daemon-reload
systemctl enable --now oligoscout-mcp
systemctl restart oligoscout-mcp
echo "== nginx: add the blocks of deploy/oligoscout-mcp.nginx.conf to the HTTPS server{} once, then:"
echo "   sudo nginx -t && sudo systemctl reload nginx"
sleep 2
curl -fsS http://127.0.0.1:8766/mcp/healthz && echo
curl -fsS http://127.0.0.1:8766/.well-known/oauth-protected-resource/mcp && echo
echo "ok - connect ChatGPT / claude.ai to https://oligoscout.org/mcp"
