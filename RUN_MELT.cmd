@echo off
rem ===========================================================================
rem  Amplikon erime egrisi + QPCR-02 (tek erime piki)
rem
rem    RUN_MELT.cmd                              qPCR tasarimi (varsayilan)
rem    RUN_MELT.cmd <tasarim_klasoru> [hedef.fa]
rem
rem  Ayarlar ortam degiskeniyle: OUT, NA, MG, DNTP, T_MIN, T_MAX, STEP,
rem  MIN_SEPARATION, TM_RESOLUTION, SEQUENCES
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; her biri asagida bash
rem   komut satirinda ACIKCA aktarilir. Yeni ayar eklerken listeye ekleyin.)
rem
rem  Grafik icin matplotlib gerekir; yoksa yalniz TSV yazilir.
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

set "DSG=%~1"
set "TGT=%~2"
if not "%DSG%"=="" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "DSG=%%j"
if not "%TGT%"=="" for /f "usebackq delims=" %%k in (`wsl.exe wslpath -a "%~f2"`) do set "TGT=%%k"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_melt.sh reference_py/cli/*.py 2>/dev/null; DESIGN='%DSG%' TARGET='%TGT%' OUT='%OUT%' NA='%NA%' MG='%MG%' DNTP='%DNTP%' T_MIN='%T_MIN%' T_MAX='%T_MAX%' STEP='%STEP%' MIN_SEPARATION='%MIN_SEPARATION%' TM_RESOLUTION='%TM_RESOLUTION%' SEQUENCES='%SEQUENCES%' bash setup/run_melt.sh"

echo.
echo  Ciktilar: amplicon_melt.tsv, amplicon_melt_curves.tsv, amplicon_melt.png
pause
