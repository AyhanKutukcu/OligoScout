@echo off
rem ===========================================================================
rem  Tek komutluk kontrol: test paketi + parametre capraz denetimi +
rem  terminal uyumsuzluk tablosunun CT-04'e etkisi
rem
rem    RUN_KONTROL.cmd
rem
rem  Ayarlar ortam degiskeniyle: OUT, TA, NA, MG, DNTP, PRIMER_NM
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; her biri asagida bash
rem   komut satirinda ACIKCA aktarilir.)
rem
rem  Hicbir adim digerini durdurmaz. Tek dosya uretir:
rem    reference_py\out\kontrol\KONTROL_RAPORU.md
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_kontrol.sh reference_py/cli/*.py 2>/dev/null; OUT='%OUT%' TA='%TA%' NA='%NA%' MG='%MG%' DNTP='%DNTP%' PRIMER_NM='%PRIMER_NM%' bash setup/run_kontrol.sh"

echo.
echo  Rapor: reference_py\out\kontrol\KONTROL_RAPORU.md
pause
