@echo off
REM =============================================================
REM  Git sync — HOME computer (D:\OneDrive\projects\df)
REM =============================================================

cd /d "D:\OneDrive\projects\df"
if errorlevel 1 (
    echo ERROR: Could not cd to D:\OneDrive\projects\df
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
