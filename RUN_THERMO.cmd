@echo off
rem ===========================================================================
rem  OligoScout - Modul 1: NN termodinamik motoru (Python referans)
rem  Cift tiklayin. Parametre defteri -> testler -> demo raporu -> panel raporu
rem  Argumanlar ortam degiskeni olarak gecer, ornegin:
rem      set TA=58 && RUN_THERMO.cmd
rem ===========================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' reference_py/RUN_THERMO.sh reference_py/cli/*.py reference_py/tests/*.py reference_py/oligoscout_thermo/*.py 2>/dev/null; TA='%TA%' NA='%NA%' MG='%MG%' DNTP='%DNTP%' PRIMER_NM='%PRIMER_NM%' TM_WINDOW='%TM_WINDOW%' PANEL='%PANEL%' bash reference_py/RUN_THERMO.sh"

echo.
echo  Ciktilar: reference_py\out\
pause
