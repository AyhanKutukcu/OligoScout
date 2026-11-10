@echo off
rem ===========================================================================
rem  Varyant katmani (VAR-01..VAR-08) — gnomAD v4
rem
rem    RUN_VARIANT.cmd                              ctDNA paneli (varsayilan)
rem    RUN_VARIANT.cmd <tarama_klasoru> <tasarim_klasoru>
rem
rem  Ornek:
rem    RUN_VARIANT.cmd reference_py\out\panel_design_qpcr_sybr ^
rem                    reference_py\out\design_qpcr_sybr
rem
rem  Ayarlar ortam degiskeniyle: SCOPE (intended|products), AF_COMMON, AF_WARN,
rem  THREE_PRIME_NT, OFFLINE, CACHE, OUT
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; her biri asagida bash
rem   komut satirinda ACIKCA aktarilir. Yeni ayar eklerken listeye ekleyin.)
rem
rem  Internet gerekir (gnomAD). Onbellek doluysa OFFLINE=1 yeter.
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

set "PNL=%~1"
set "DSG=%~2"
if not "%PNL%"=="" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "PNL=%%j"
if not "%DSG%"=="" for /f "usebackq delims=" %%k in (`wsl.exe wslpath -a "%~f2"`) do set "DSG=%%k"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_variant.sh reference_py/cli/*.py 2>/dev/null; PANEL='%PNL%' DESIGN='%DSG%' SCOPE='%SCOPE%' AF_COMMON='%AF_COMMON%' AF_WARN='%AF_WARN%' THREE_PRIME_NT='%THREE_PRIME_NT%' OFFLINE='%OFFLINE%' CACHE='%CACHE%' OUT='%OUT%' bash setup/run_variant.sh"

echo.
echo  Ciktilar: variant_sites.tsv, variant_qc.tsv, variant_qc_records.jsonl
pause
