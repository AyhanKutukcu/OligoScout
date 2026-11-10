@echo off
rem  RT-03 izoform kapsami (GENCODE ekson anotasyonu; transkriptom FM-index YOK)
rem    set PAIRS=C:\yol\pairs.tsv
rem    RUN_ISOFORM.cmd
setlocal
cd /d "%~dp0"
wsl.exe -e true >nul 2>&1
if errorlevel 1 ( echo [HATA] WSL yok. RUN_SETUP.cmd & pause & exit /b 1 )
for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
set "PRS=%PAIRS%"
if not "%PRS%"=="" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%PRS%"`) do set "PRS=%%j"
wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_isoform.sh reference_py/cli/*.py 2>/dev/null; PAIRS='%PRS%' EXONS='%EXONS%' OUT='%OUT%' bash setup/run_isoform.sh"
echo. & echo  Cikti: reference_py\out\isoform\isoform_qc.tsv
pause
