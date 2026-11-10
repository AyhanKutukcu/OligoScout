@echo off
rem ===========================================================================
rem  OligoScout - TEK GIRIS NOKTASI
rem  Bu dosyaya cift tiklayin. Her sey WSL icinde sirayla calisir.
rem  Argumanlar dogrudan setup\00_bootstrap.sh scriptine gecer, ornegin:
rem      RUN_SETUP.cmd --from 5
rem      RUN_SETUP.cmd --only 8
rem ===========================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

where wsl.exe >nul 2>&1
if errorlevel 1 (
    echo [HATA] wsl.exe bulunamadi. Yonetici PowerShell'de: wsl --install -d Ubuntu
    pause & exit /b 1
)

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok.
    echo        Yonetici PowerShell'de:  wsl --install -d Ubuntu
    echo        Kurulumdan sonra bir kez Ubuntu'yu acip kullanici adi/parola belirleyin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
if "%WSLDIR%"=="" (
    echo [HATA] Klasor yolu WSL'e cevrilemedi.
    pause & exit /b 1
)

echo.
echo  Calisma dizini (WSL) : %WSLDIR%
echo  Baslatiliyor         : setup/00_bootstrap.sh %*
echo.

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/*.sh setup/*.py 2>/dev/null; bash setup/00_bootstrap.sh %*"

echo.
echo  Bitti. Loglar: logs\
pause
