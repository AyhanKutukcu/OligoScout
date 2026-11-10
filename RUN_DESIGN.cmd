@echo off
rem ===========================================================================
rem  Primer tasarimi (ozgulluk taramasi ONCESI)
rem
rem    RUN_DESIGN.cmd <hedef.fa>
rem    RUN_DESIGN.cmd <hedef.fa> <profil>
rem    RUN_DESIGN.cmd ACGTACGT...
rem
rem  Ayarlar ortam degiskeniyle: PROFILE, PRODUCT, TM, LENGTH, GC, TA, MG,
rem                              TOP, PROBE, PROBE_LENGTH, TM_DELTA, PROBES,
rem                              REQUIRE_TH03 (1: TH-03 FAIL primerlerini ele)
rem  (WSL bu degiskenleri Windows'tan devralmaz; her biri asagida bash
rem   komut satirinda ACIKCA aktarilir. Yeni ayar eklerken listeye ekleyin.)
rem
rem  Hidroliz probu (PRB-01..PRB-05):
rem    RUN_DESIGN.cmd hedef.fa ddpcr      profil probu kendiliginden secer
rem    RUN_DESIGN.cmd hedef.fa qpcr_sybr prob   ucuncu arguman "prob" ise
rem                                             prob adimi zorlanir
rem    set PRODUCT=70:150 && RUN_DESIGN.cmd hedef.fa
rem
rem  Uygulama profili (kisitlari kendisi belirler) - ikinci arguman:
rem    RUN_DESIGN.cmd hedef.fa ctdna
rem  (PowerShell'de:  .\RUN_DESIGN.cmd hedef.fa ctdna )
rem  Profiller: qpcr_sybr rt_qpcr ctdna ddpcr sanger ngs_amplicon
rem             standard_pcr bisulfite
rem
rem  Sonra ozgulluk icin:  set PANEL=... && RUN_PANEL.cmd
rem ===========================================================================
setlocal
cd /d "%~dp0"

if "%~1"=="" (
    echo Kullanim: RUN_DESIGN.cmd ^<hedef.fa ^| dizi^>
    pause & exit /b 2
)

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

rem Ikinci arguman verilmisse profil odur; yoksa ortam degiskeni gecerli.
if not "%~2"=="" set "PROFILE=%~2"
rem Ucuncu arguman "prob" ise prob modulu acilir.
if /i "%~3"=="prob" set "PROBE=1"

set "TGT=%~1"
if exist "%~1" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "TGT=%%j"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_design.sh reference_py/cli/*.py 2>/dev/null; TARGET='%TGT%' PROFILE='%PROFILE%' PRODUCT='%PRODUCT%' TM='%TM%' LENGTH='%LENGTH%' GC='%GC%' TA='%TA%' MG='%MG%' TOP='%TOP%' PROBE='%PROBE%' PROBE_LENGTH='%PROBE_LENGTH%' TM_DELTA='%TM_DELTA%' PROBES='%PROBES%' REQUIRE_TH03='%REQUIRE_TH03%' bash setup/run_design.sh"

echo.
if "%PROFILE%"=="" (
    echo  Ciktilar: reference_py\out\design\
) else (
    echo  Ciktilar: reference_py\out\design_%PROFILE%\
)
pause
