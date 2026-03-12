@echo off
setlocal

REM --- Paths ---
set "R_SCRIPT=C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
set "SCRIPT=D:\OneDrive - NRCan RNCan\gs\doug-fir\snp_qc_pipeline_v2.R"
set "OUT_DIR=D:\OneDrive - NRCan RNCan\gs\doug-fir\data\qc"

REM --- Timestamp via PowerShell (safe on Windows) ---
for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"') do set "STAMP=%%i"

set "LOG_FILE=%OUT_DIR%\qc_report_%STAMP%.txt"

if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"

echo ======================================================
echo Running SNP QC pipeline...
echo Script: %SCRIPT%
echo Log   : %LOG_FILE%
echo ======================================================

"%R_SCRIPT%" "%SCRIPT%" > "%LOG_FILE%" 2>&1

set EXITCODE=%ERRORLEVEL%

echo.
echo ======================================================
if %EXITCODE%==0 (
    echo Pipeline finished successfully.
) else (
    echo Pipeline failed with exit code %EXITCODE%.
)
echo Log file: %LOG_FILE%
echo ======================================================

pause
exit /b %EXITCODE%