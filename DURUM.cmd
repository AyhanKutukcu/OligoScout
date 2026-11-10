@echo off
rem ===========================================================================
rem  PPFM indeks uretiminin canli durumu.
rem  Cift tiklayin: 10 saniyede bir yenilenir. Cikmak icin Ctrl+C.
rem  Tek seferlik rapor icin:  DURUM.cmd tek
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

set "ARG=-w"
if /i "%~1"=="tek" set "ARG="

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/ppfm_status.sh 2>/dev/null; bash setup/ppfm_status.sh %ARG%"
