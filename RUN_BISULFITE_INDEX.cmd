@echo off
rem  Bisulfit-donusturulmus genom indeksleri (BIS-02/03/04) - AGIR, saatler surer
rem    RUN_BISULFITE_INDEX.cmd
setlocal
cd /d "%~dp0"
wsl.exe -e true >nul 2>&1
if errorlevel 1 ( echo [HATA] WSL yok. RUN_SETUP.cmd & pause & exit /b 1 )
for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_bisulfite_index.sh 2>/dev/null; PRESET='%PRESET%' SA_RATE='%SA_RATE%' CHROMS='%CHROMS%' JOBS='%JOBS%' bash setup/run_bisulfite_index.sh"
echo. & echo  Cikti: OligoScout\data\indexes\bisulfite\{CT,GA}
pause
