@echo off
rem ===========================================================================
rem  OligoScout panel boru hatti
rem  Derle -> tum genomda panel aramasi -> termodinamik QC birlestirmesi
rem
rem    RUN_PANEL.cmd                     depodaki chr22 32-cift paneli
rem    RUN_PANEL.cmd <panel.tsv>         kendi paneliniz (Windows yolu olur)
rem    RUN_PANEL.cmd <panel.tsv> 166     ctDNA: urunler cfDNA fragmaniyla
rem                                      agirliklanir (MPX-02f)
rem    RUN_PANEL.cmd <panel.tsv> 166 hizli   taramayi atla, yalniz QC'yi
rem                                      yeniden hesapla (saniyeler)
rem
rem  Panel basligi: pair_id + (primer1/primer2 ya da forward_primer/reverse_primer)
rem  Ayarlar ortam degiskeniyle: TA, NA, MG, DNTP, PRIMER_NM, TM_WINDOW,
rem  ANCHOR, MAXMM, MINAMP, MAXAMP, CACHE, FRAGMENT_BP, SKIP_SCAN
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; her biri asagida
rem   bash komut satirinda ACIKCA aktarilir. Listeye eklemeyi unutmayin.)
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

rem Panel: once komut satiri argumani, yoksa PANEL ortam degiskeni.
set "PANELARG=%~1"
if not "%PANELARG%"=="" (
    if exist "%~1" (
        for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "PANEL=%%j"
    ) else (
        echo [HATA] Panel dosyasi bulunamadi: %~1
        pause & exit /b 2
    )
)

rem Ikinci arguman: cfDNA fragman uzunlugu. Ucuncu: "hizli" -> taramayi atla.
if not "%~2"=="" set "FRAGMENT_BP=%~2"
if /i "%~3"=="hizli" set "SKIP_SCAN=1"

if not "%PANEL%"=="" echo  Panel: %PANEL%
if not "%FRAGMENT_BP%"=="" echo  cfDNA fragman: %FRAGMENT_BP% bp
if "%SKIP_SCAN%"=="1" echo  Tarama atlaniyor; mevcut baglanma kayitlari kullanilacak.

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_panel.sh reference_py/cli/*.py 2>/dev/null; PANEL='%PANEL%' TA='%TA%' NA='%NA%' MG='%MG%' DNTP='%DNTP%' PRIMER_NM='%PRIMER_NM%' TM_WINDOW='%TM_WINDOW%' ANCHOR='%ANCHOR%' MAXMM='%MAXMM%' MINAMP='%MINAMP%' MAXAMP='%MAXAMP%' CACHE='%CACHE%' FRAGMENT_BP='%FRAGMENT_BP%' SKIP_SCAN='%SKIP_SCAN%' bash setup/run_panel.sh"

echo.
echo  Ciktilar yukarida listelenen klasorde.
pause
