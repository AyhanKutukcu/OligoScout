@echo off
rem ===========================================================================
rem  Dizi/kompozisyon/cift/amplikon kontrolleri (SEQ, CMP, SS-02, PAIR, AMP)
rem
rem    RUN_SEQUENCE_QC.cmd <tasarim_klasoru> [hedef.fa] [profil]
rem    RUN_SEQUENCE_QC.cmd                     (varsayilan: ctdna tasarimi)
rem
rem  Ayarlar ortam degiskeniyle: OUT, NA, MG, DNTP, PRIMER_NM, TA, PRODUCT_BP
rem  Ciktilar: seqpair_qc.tsv, seqpair_panel.tsv, seqpair_qc_records.jsonl
rem ===========================================================================
setlocal
cd /d "%~dp0"
wsl.exe -e true >nul 2>&1
if errorlevel 1 (
    echo [HATA] Calisan bir WSL dagitimi yok. Once RUN_SETUP.cmd calistirin.
    pause & exit /b 1
)
for /f "usebackq delims=" %%i in (`wsl.exe wslpath -a "%~dp0."`) do set "WSLDIR=%%i"
set "DSG=%~1"
set "TGT=%~2"
set "PRF=%~3"
if not "%DSG%"=="" for /f "usebackq delims=" %%j in (`wsl.exe wslpath -a "%~f1"`) do set "DSG=%%j"
if not "%TGT%"=="" for /f "usebackq delims=" %%k in (`wsl.exe wslpath -a "%~f2"`) do set "TGT=%%k"
wsl.exe -e bash -lc "cd '%WSLDIR%' && sed -i 's/\r$//' setup/run_sequence_qc.sh reference_py/cli/*.py 2>/dev/null; DESIGN='%DSG%' TARGET='%TGT%' PROFILE='%PRF%' GENE='%GENE%' PROBES='%PROBES%' OUT='%OUT%' NA='%NA%' MG='%MG%' DNTP='%DNTP%' PRIMER_NM='%PRIMER_NM%' TA='%TA%' PRODUCT_BP='%PRODUCT_BP%' bash setup/run_sequence_qc.sh"
echo.
echo  Ciktilar: seqpair_qc.tsv, seqpair_panel.tsv
pause
