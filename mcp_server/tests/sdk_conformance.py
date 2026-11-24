import asyncio, json, os, subprocess, sys, time
import mcp
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
# Conformance check with the OFFICIAL MCP Python SDK (pip install mcp).
#   python3 mcp_server/tests/sdk_conformance.py
# Exercises stdio and Streamable HTTP: 13 checks each (26 total).
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("OLIGOSCOUT_ROOT") or os.path.dirname(os.path.dirname(HERE))
TPL = os.path.join(HERE, "tpl.fa")
CHECKS = []
def ok(name, cond, info=""):
    CHECKS.append((name, bool(cond), info)); print(("PASS " if cond else "FAIL ") + name, info)

async def exercise(session, label):
    init = session._server_info if hasattr(session, "_server_info") else None
    tools = (await session.list_tools()).tools
    ok(f"{label}: tools/list parsed by SDK", len(tools) == 15, f"{len(tools)} tools")
    ok(f"{label}: every tool has inputSchema+annotations", all(t.input_schema and t.annotations for t in tools))
    r = await session.call_tool("oligoscout_thermo", {"sequences": ["TCCTCAAAGGCACCTACTG"]})
    sc = r.structured_content or {}
    ok(f"{label}: thermo Tm == engine record 60.1601", abs(sc["results"][0]["tm_c"] - 60.1601) < 1e-4, str(sc["results"][0]["tm_c"]))
    ok(f"{label}: provenance carries REP-03 + digest", "REP-03" in sc["provenance"]["calibration"] and len(sc["provenance"]["param_set_sha256"]) == 64)
    r = await session.call_tool("oligoscout_dimer", {"sequence_a": "CTCCTCAAAGGCACCTACT", "sequence_b": "GGGAACACTAGCATTCACAC"})
    sys.path.insert(0, os.path.join(ROOT, "reference_py"))
    from oligoscout_thermo import align
    want = align.best_dimer("CTCCTCAAAGGCACCTACT", "GGGAACACTAGCATTCACAC").dG37_kcal
    ok(f"{label}: hetero-dimer ΔG == align.best_dimer", abs(r.structured_content["dimers"]["any"]["dG37_kcal_mol"] - want) < 1e-4, f"{want:.4f}")
    r = await session.call_tool("oligoscout_check_pair", {"forward": "ACGT", "reverse": "ACGTACGTACGTAC"})
    ok(f"{label}: schema violation -> isError tool result", r.is_error and "Girdi" in r.content[0].text, r.content[0].text[:60])
    r = await session.call_tool("oligoscout_check_pair", {"forward": "CTCCTCAAAGGCACCTACT", "reverse": "GGGAACACTAGCATTCACAC",
        "template_sequence": open(TPL).read(), "profile": "ctdna"})
    sc = r.structured_content
    ok(f"{label}: pair QC on pasted FASTA", not r.is_error and sc["pair"]["product_bp"] == 71 and "VAR-01..08" in sc["not_evaluated"], f"worst {sc.get('worst_severity')} {sc.get('counts')}")
    r = await session.call_tool("oligoscout_design_primers", {"template_sequence": open(TPL).read(), "profile": "ctdna", "top": 3, "specificity": False})
    sc = r.structured_content
    ok(f"{label}: design on FASTA -> top pair == design_primers.py rank 1", sc["pairs"][0]["forward"] == "TCCTCAAAGGCACCTACTG" and sc["pairs"][0]["reverse"] == "AGCATTCACACCGTAGCA",
       f"{sc['pairs'][0]['forward']}/{sc['pairs'][0]['reverse']} · not_evaluated={list(sc['not_evaluated'])}")
    r = await session.call_tool("oligoscout_explain_rule", {"rule_id": "VAR-02"})
    ok(f"{label}: explain_rule", r.structured_content["rules"][0]["id"] == "VAR-02")
    res = (await session.list_resources()).resources
    ok(f"{label}: resources/list", {str(x.uri) for x in res} >= {"oligoscout://guide", "oligoscout://catalog"})
    txt = (await session.read_resource("oligoscout://rule/TH-02")).contents[0].text
    ok(f"{label}: resource template rule/TH-02", json.loads(txt)["id"] == "TH-02")
    pr = await session.get_prompt("design_assay", {"target": "EGFR exon 19", "application": "ddPCR"})
    ok(f"{label}: prompts/get", "EGFR" in pr.messages[0].content.text)

async def main():
    params = StdioServerParameters(command=sys.executable, args=[f"{ROOT}/mcp_server/server.py"], env={**os.environ, "OLIGOSCOUT_GNOMAD_API": "0"})
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            info = await s.initialize()
            ok("stdio: initialize", info.server_info.name == "oligoscout", f"protocol {info.protocol_version}")
            await exercise(s, "stdio")
    proc = subprocess.Popen([sys.executable, f"{ROOT}/mcp_server/server.py", "--http", "--port", "8799"], env={**os.environ, "OLIGOSCOUT_GNOMAD_API": "0"}, stderr=subprocess.DEVNULL)
    time.sleep(2.5)
    try:
        try:
            from mcp.client.streamable_http import streamablehttp_client as shc
        except ImportError:
            from mcp.client.streamable_http import streamable_http_client as shc
        async with shc("http://127.0.0.1:8799/mcp") as streams:
            r, w = streams[0], streams[1]
            async with ClientSession(r, w) as s:
                info = await s.initialize()
                ok("http: initialize", info.server_info.name == "oligoscout", f"protocol {info.protocol_version}")
                await exercise(s, "http")
    finally:
        proc.terminate()
    print(f"\n{sum(c[1] for c in CHECKS)}/{len(CHECKS)} checks passed")
    json.dump([{"check": c[0], "pass": c[1], "info": c[2]} for c in CHECKS], open(os.path.join(HERE, "sdk_results.json"), "w"), ensure_ascii=False, indent=1)
asyncio.run(main())
