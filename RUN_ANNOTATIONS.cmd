@echo off
rem ===========================================================================
rem  Anotasyon on-islemesi (SPC-06 psodogen, SPC-07 segmental duplikasyon)
rem
rem    RUN_ANNOTATIONS.cmd
rem
rem  OligoScout\data\annotations icindeki ham GTF (GENCODE) ve
rem  genomicSuperDups dosyalarini kompakt aralik tablolarina indirger.
rem  Bir kez calistirin; panel_qc bu tablolari otomatik kullanir.
rem ===========================================================================
setlocal
cd /d "%~dp0"
wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)
for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' reference_py/cli/build_annotations.py 2>/dev/null; cd reference_py && python3 cli/build_annotations.py"
echo.
echo  Ciktilar: OligoScout\data\annotations\pseudogene_intervals.tsv, segdup_intervals.tsv
pause
