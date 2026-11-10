@echo off
rem ===========================================================================
rem  Iki profilin karsilastirma sekli
rem
rem    RUN_COMPARE.cmd                       varsayilan iki profil
rem    RUN_COMPARE.cmd <klasorA> <klasorB>   kendi tarama klasorleriniz
rem
rem  Girdi: RUN_PANEL.cmd'nin yazdigi tarama klasorleri (panel_qc.tsv iceren).
rem  Cikti: reference_py\out\profile_comparison.png ve .tsv
rem
rem  Ayarlar ortam degiskeniyle: NAME_A, NAME_B, LABEL_A, LABEL_B, OUT
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; her biri asagida bash
rem   komut satirinda ACIKCA aktarilir. Yeni ayar eklerken listeye ekleyin.)
rem
rem  Bu adim WSL icinde matplotlib + numpy ister; eksikse kurulum satirini yazar.
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

set "PA=%~1"
set "PB=%~2"
if not "%PA%"=="" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "PA=%%j"
if not "%PB%"=="" for /f "usebackq delims=" %%k in (`wsl.exe wslpath -a "%~f2"`) do set "PB=%%k"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_compare.sh reference_py/cli/*.py 2>/dev/null; PANEL_A='%PA%' PANEL_B='%PB%' NAME_A='%NAME_A%' NAME_B='%NAME_B%' LABEL_A='%LABEL_A%' LABEL_B='%LABEL_B%' OUT='%OUT%' bash setup/run_compare.sh"

echo.
echo  Ciktilar: reference_py\out\profile_comparison.png ve .tsv
pause
