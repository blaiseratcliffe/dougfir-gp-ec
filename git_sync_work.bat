@echo off
REM =============================================================
REM  Git sync — WORK computer (D:\OneDrive - NRCan RNCan\gs\doug-fir)
REM =============================================================

cd /d "D:\OneDrive - NRCan RNCan\gs\doug-fir"
if errorlevel 1 (
    echo ERROR: Could not cd to D:\OneDrive - NRCan RNCan\gs\doug-fir
    pause
    exit /b 1
)

echo.
echo === Pull latest changes ===
git pull
echo.

echo === Current status ===
git status --short
echo.

REM --- If there are no local changes (tracked or untracked), stop here ---
for /f %%i in ('git status --porcelain') do goto HAS_CHANGES
echo No local changes to commit.
pause
exit /b 0
:HAS_CHANGES

echo === Files to be committed: ===
git add -A --dry-run
echo.
set /p MSG="Commit message (or Ctrl+C to cancel): "
git add -A
git commit -m "%MSG%"
git push

echo.
echo === Done. ===
pause
