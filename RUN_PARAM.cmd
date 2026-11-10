@echo off
rem ===========================================================================
rem  Parametre aktarim is akisi
rem
rem    RUN_PARAM.cmd              params_in\ icine bos sablonlar uretir
rem    RUN_PARAM.cmd check        doldurdugunuz dosyalari denetler
rem    RUN_PARAM.cmd emit         params.py'ye yapistirilacak bloklari yazar
rem
rem  Cift aktarim icin ikinci kisinin dosyasini params_in\<tablo>.2.tsv
rem  olarak kaydedin; check ikisini otomatik karsilastirir.
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

set "SUB=%~1"
if "%SUB%"=="" set "SUB=template"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/params.sh reference_py/cli/*.py 2>/dev/null; bash setup/params.sh %SUB%"

echo.
pause
