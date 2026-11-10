@echo off
rem ===========================================================================
rem  Tekrarlanabilir uctan uca regresyon
rem
rem    RUN_REGRESYON.cmd
rem    RUN_REGRESYON.cmd <selected_panel.tsv>    (TH-02/TH-04 adimini da kosar)
rem
rem  Ayarlar ortam degiskeniyle: OUT, TIMEOUT, HIZLI=1
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; asagida ACIKCA aktarilir.)
rem
rem  Test paketini ve 16 CLI'nin tamamini ayni girdilerle siraya dizer,
rem  ciktilarin sema basliklarini denetler, tek rapor yazar:
rem    reference_py\out\regresyon\REGRESYON_RAPORU.md
rem
rem  Girdisi olmayan adimlar FAIL degil SKIP olur (ornegin C++ dokumu yoksa).
rem  Herhangi bir adim FAIL ise cikis kodu 1.
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

set "LST=%~1"
if not "%LST%"=="" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "LST=%%j"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_regresyon.sh reference_py/cli/*.py 2>/dev/null; OUT='%OUT%' TIMEOUT='%TIMEOUT%' HIZLI='%HIZLI%' LATEST='%LST%' bash setup/run_regresyon.sh"

echo.
echo  Rapor: reference_py\out\regresyon\REGRESYON_RAPORU.md
pause
