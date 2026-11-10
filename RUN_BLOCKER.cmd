@echo off
rem  CT-05 (WT blocker / ddPCR prob ayrimi) + prob x primer dimerleri
rem    RUN_BLOCKER.cmd  (diziler asagida; kendi reaktiflerinizle degistirin)
rem  LNA stabilizasyonu MODELLENMEZ; Tm'ler DNA iskeleti alt siniri.
setlocal
cd /d "%~dp0"
wsl.exe -e true >nul 2>&1
if errorlevel 1 ( echo [HATA] WSL yok. RUN_SETUP.cmd & pause & exit /b 1 )
for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
wsl.exe -e bash -lc "cd '%WSLDIR%/reference_py' && sed -i 's/\r$//' cli/*.py 2>/dev/null; python3 cli/blocker_qc.py --mutant 'CCC+CT+C+T+C+TGAC' --wt 'C+CT+C+T+ATG+A+CC' --forward 'AGTCCTTTCTGTGTCTTCCCA' --reverse 'TCCAGCATCTCCAGCAGC' --blocker 'AGCATCTCCAGCAGCAGG+T+CA+T+AGAGGG/3Phos/' --out out/blocker"
echo. & echo  Cikti: reference_py\out\blocker\ct05_qc.tsv
pause
