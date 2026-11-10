@echo off
rem  BIS-04 bisulfit iplik ozgullugu - donusturulmus genom taramasi
rem    RUN_BISULFITE_PANEL.cmd
setlocal
cd /d "%~dp0"
wsl.exe -e true >nul 2>&1
if errorlevel 1 ( echo [HATA] WSL yok. RUN_SETUP.cmd & pause & exit /b 1 )
for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_bisulfite_panel.sh 2>/dev/null; PANEL='%PANEL%' INTENDED='%INTENDED%' bash setup/run_bisulfite_panel.sh"
echo. & echo  Cikti: reference_py\out\bisulfite_panel\bisulfite_qc.tsv
pause
