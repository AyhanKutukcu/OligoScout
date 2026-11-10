@echo off
rem ===========================================================================
rem  Python <-> C++ diferansiyel testi (native_thermo cekirdegi)
rem
rem    RUN_DIFF_CPP.cmd
rem
rem  Ayarlar ortam degiskeniyle: PRESET, SEQUENCES, OUT, NA, MG, DNTP,
rem  PRIMER_NM, TA, SALT_MODEL
rem  (WSL ortam degiskenlerini Windows'tan devralmaz; her biri asagida bash
rem   komut satirinda ACIKCA aktarilir. Yeni ayar eklerken listeye ekleyin.)
rem
rem  Once `oligoscout_thermo_dump` hedefini derler. Derleme yapilandirmasi
rem  yoksa once RUN_SETUP.cmd --only 5 calistirin.
rem
rem  Cekirdek sutunlarda fark cikarsa cikis kodu 1 olur. Dimer sutunlarindaki
rem  fark beklenir ve ayri sinif olarak raporlanir.
rem ===========================================================================
setlocal
cd /d "%~dp0"

wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"

wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_diff_cpp.sh reference_py/cli/*.py 2>/dev/null; PRESET='%PRESET%' SEQUENCES='%SEQUENCES%' OUT='%OUT%' NA='%NA%' MG='%MG%' DNTP='%DNTP%' PRIMER_NM='%PRIMER_NM%' TA='%TA%' SALT_MODEL='%SALT_MODEL%' bash setup/run_diff_cpp.sh"

echo.
echo  Ciktilar: reference_py\out\diff_cpp\diff_cpp.tsv, diff_cpp_summary.md
pause
