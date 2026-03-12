@echo off
REM =============================================================
REM  Initialize git in the existing working directory
REM  
REM  BEFORE RUNNING:
REM    1. Place .gitignore and README.md in D:\OneDrive - NRCan RNCan\gs\doug-fir\
REM    2. Edit the GitHub remote URL below
REM    3. Double-click this script (or run from cmd in that directory)
REM =============================================================

echo.
echo === Douglas-fir GS repository setup ===
echo.

cd /d "D:\OneDrive - NRCan RNCan\gs\doug-fir"
if errorlevel 1 (
    echo ERROR: Could not cd to D:\OneDrive - NRCan RNCan\gs\doug-fir
    pause
    exit /b 1
)

REM --- Initialize ---
git init
git branch -M main

REM --- Dry run: see what git will track ---
echo.
echo === Files that will be tracked (dry run): ===
git add -A --dry-run
echo.
echo Review the list above. If anything unexpected appears, press Ctrl+C now.
pause

REM --- Stage everything the whitelist allows ---
git add -A

REM --- Initial commit ---
git commit -m "Initial commit: QC pipeline, site analysis, and project scripts"

REM --- Connect to GitHub remote ---
REM Uncomment and edit the two lines below, then re-run the script:
REM git remote add origin https://github.com/YOUR_USERNAME/dfir-gs.git
REM git push -u origin main

echo.
echo === Done. ===
echo   - Run 'git log --oneline' to verify the commit
echo   - Run 'git status' to confirm nothing unexpected is untracked
echo.
pause
