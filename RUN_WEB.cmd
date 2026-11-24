@echo off
setlocal
pushd "%~dp0"
wsl.exe -d Ubuntu --exec bash -lc "python3 web/build_web_data.py && python3 web/build_site.py && exec python3 web/local_server.py"
set "OLIGOSCOUT_WEB_EXIT=%errorlevel%"
popd
if not "%OLIGOSCOUT_WEB_EXIT%"=="0" pause
exit /b %OLIGOSCOUT_WEB_EXIT%
