@echo off
setlocal EnableExtensions

cd /d "%~dp0"

title Candidate Gene SNP Annotation Pipeline

echo.
echo ============================================================
echo        CANDIDATE GENE SNP ANNOTATION PIPELINE
echo ============================================================
echo Project: %CD%
echo.

REM ------------------------------------------------------------
REM Fixed project inputs
REM ------------------------------------------------------------
set "INPUT_FILE=SNP_list_common.csv"
set "REFERENCE_DIR=references"
set "REFERENCE_FASTA=%REFERENCE_DIR%\reference.fa"
set "ANNOTATION_GFF=%REFERENCE_DIR%\annotation.gff3"
set "OUTPUT_DIR=results"

REM ------------------------------------------------------------
REM Check Python and launcher
REM ------------------------------------------------------------
python --version
if errorlevel 1 (
    echo ERROR: Python was not found in PATH.
    pause
    exit /b 1
)

if not exist "12_run_pipeline.py" (
    echo ERROR: 12_run_pipeline.py not found.
    pause
    exit /b 1
)

if not exist "modules\pipeline_common.py" (
    echo ERROR: modules\pipeline_common.py not found.
    pause
    exit /b 1
)

REM ------------------------------------------------------------
REM Check user inputs
REM ------------------------------------------------------------
if not exist "%INPUT_FILE%" (
    echo ERROR: %INPUT_FILE% not found.
    pause
    exit /b 1
)

if not exist "%REFERENCE_FASTA%" (
    echo ERROR: %REFERENCE_FASTA% not found.
    echo Put the reference FASTA inside the references folder.
    pause
    exit /b 1
)

if not exist "%ANNOTATION_GFF%" (
    echo ERROR: %ANNOTATION_GFF% not found.
    echo Put the GFF/GTF annotation inside the references folder.
    pause
    exit /b 1
)

if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

echo.
echo ============================================================
echo INPUTS
echo ============================================================
echo SNP file   : %INPUT_FILE%
echo Reference  : %REFERENCE_FASTA%
echo Annotation : %ANNOTATION_GFF%
echo Output      : %OUTPUT_DIR%\
echo ============================================================
echo.

REM ------------------------------------------------------------
REM Test launcher
REM ------------------------------------------------------------
python 12_run_pipeline.py --help
if errorlevel 1 (
    echo.
    echo ERROR: Pipeline launcher test failed.
    pause
    exit /b 1
)

echo.
echo Launcher test passed.
echo.

REM ------------------------------------------------------------
REM Run complete pipeline
REM ------------------------------------------------------------
python 12_run_pipeline.py ^
    --input "%INPUT_FILE%" ^
    --reference "%REFERENCE_FASTA%" ^
    --gff "%ANNOTATION_GFF%" ^
    --out "%OUTPUT_DIR%" ^
    --species "triticum_aestivum" ^
    --window 100000 ^
    --min-identity 90 ^
    --min-coverage 70 ^
    --flank 250 ^
    --online-ncbi

if errorlevel 1 (
    echo.
    echo ============================================================
    echo                 PIPELINE FAILED
    echo ============================================================
    echo.
    pause
    exit /b 1
)

echo.
echo ============================================================
echo              PIPELINE COMPLETED
 echo ============================================================
echo Final workbook:
echo %OUTPUT_DIR%\candidate_gene_master.xlsx
echo.
pause
