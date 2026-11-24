@echo off
rem OligoScout MCP sunucusu (Streamable HTTP + REST + OpenAPI) - http://127.0.0.1:8766/mcp
rem Claude Desktop / Cursor gibi yerel istemciler bu dosyaya ihtiyac duymaz (stdio, bkz. mcp_server\README.md).
setlocal
pushd "%~dp0"
wsl.exe -d Ubuntu --exec bash -lc "exec python3 mcp_server/server.py --http --host 127.0.0.1 --port 8766"
set "OLIGOSCOUT_MCP_EXIT=%errorlevel%"
popd
if not "%OLIGOSCOUT_MCP_EXIT%"=="0" pause
exit /b %OLIGOSCOUT_MCP_EXIT%
