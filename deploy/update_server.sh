#!/usr/bin/env bash
# One-step update of oligoscout.org (web engine + MCP server) after a GitHub push.
#
# Run ON THE SERVER from any directory:
#     sudo bash /opt/oligoscout/deploy/update_server.sh
#
# Steps: git pull (as the service user) -> C++ tools -> tests -> SPC-08 extra contigs
# (once) -> web data + single-file site -> systemd units -> restart -> health checks.
# It never touches credentials; the optional MCP token stays in /etc/oligoscout/mcp.env.
# nginx is checked, not edited: the MCP blocks are added once by hand (see the end).
set -euo pipefail
APP=${APP:-/opt/oligoscout}
SVC_USER=${SVC_USER:-oligoscout}
as_app() { sudo -u "$SVC_USER" -H bash -c "cd '$APP' && $*"; }

[[ $EUID -eq 0 ]] || { echo "run with sudo"; exit 1; }
[[ -d "$APP/.git" ]] || { echo "no checkout at $APP"; exit 1; }

echo "== 1/7 git pull"
as_app "git fetch --tags origin && git checkout main && git pull --ff-only origin main"
as_app "git log --oneline -1"

echo "== 2/7 C++ tools (panel, contig scan, virtual PCR scan)"
as_app "cd OligoScout && cmake --preset release >/dev/null && \
        cmake --build --preset release --target oligoscout_panel oligoscout_contig_scan \
              oligoscout_vpcr_scan test_single_primer_products test_terminal_rule_search -j\$(nproc) && \
        ./build/release/test_single_primer_products && ./build/release/test_terminal_rule_search"

echo "== 3/7 Python and MCP tests"
# generated, git-ignored catalogue outputs (reference_py/out/durum) follow the YAML catalogue
as_app "python3 reference_py/cli/catalog_yaml.py export >/dev/null"
as_app "cd reference_py && python3 -B -m unittest discover -s tests -q"
as_app "python3 -B -m unittest discover -s mcp_server/tests -q"

echo "== 4/7 SPC-08 extra contigs (built once from the RefSeq assembly)"
if [[ ! -s "$APP/OligoScout/data/annotations/grch38_extra_contigs.fa" ]]; then
    if [[ -s "$APP/downloads/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" ]]; then
        as_app "python3 reference_py/cli/build_extra_contigs.py"
    else
        echo "   skipped: downloads/GCF_000001405.40_GRCh38.p14_genomic.fna.gz not found (SPC-08 falls back to chrM)"
    fi
fi

echo "== 5/7 web data + site"
as_app "python3 web/build_web_data.py && python3 web/build_site.py"

echo "== 6/7 systemd units + restart"
install -d -m 750 -o "$SVC_USER" -g "$SVC_USER" /etc/oligoscout
[[ -f /etc/oligoscout/mcp.env ]] || install -m 600 -o "$SVC_USER" -g "$SVC_USER" /dev/null /etc/oligoscout/mcp.env
cp "$APP/deploy/oligoscout.service" /etc/systemd/system/oligoscout.service
cp "$APP/deploy/oligoscout-mcp.service" /etc/systemd/system/oligoscout-mcp.service
systemctl daemon-reload
systemctl enable oligoscout oligoscout-mcp >/dev/null
systemctl restart oligoscout oligoscout-mcp
sleep 3

echo "== 7/7 health checks"
ok=1
curl -fsS -o /dev/null http://127.0.0.1:8765/ && echo "   web     : ok (127.0.0.1:8765)" || { echo "   web     : FAILED"; ok=0; }
curl -fsS http://127.0.0.1:8766/mcp/healthz >/dev/null && echo "   mcp     : ok (127.0.0.1:8766)" || { echo "   mcp     : FAILED"; ok=0; }
curl -fsS http://127.0.0.1:8766/.well-known/oauth-protected-resource/mcp >/dev/null \
    && echo "   oauth   : ok (protected-resource metadata)" || { echo "   oauth   : FAILED"; ok=0; }
if nginx -T 2>/dev/null | grep -q "location = /mcp"; then
    nginx -t >/dev/null 2>&1 && systemctl reload nginx && echo "   nginx   : MCP blocks present, reloaded"
    curl -fsS -o /dev/null https://oligoscout.org/mcp/healthz && echo "   public  : https://oligoscout.org/mcp ok" \
        || echo "   public  : https://oligoscout.org/mcp not reachable yet (DNS/certificate?)"
else
    echo "   nginx   : MCP blocks NOT found. Add the blocks of $APP/deploy/oligoscout-mcp.nginx.conf"
    echo "             inside the HTTPS server { } of oligoscout.org, then: sudo nginx -t && sudo systemctl reload nginx"
fi
journalctl -u oligoscout -u oligoscout-mcp -n 5 --no-pager || true
[[ $ok -eq 1 ]] && echo "done - site: https://oligoscout.org  MCP: https://oligoscout.org/mcp" || exit 1
