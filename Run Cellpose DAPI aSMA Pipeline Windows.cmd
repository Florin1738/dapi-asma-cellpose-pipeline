@echo off
setlocal

set "PROJECT_DIR=%~dp0"
rem %~dp0 always ends in a backslash; a trailing "\" before the closing quote
rem escapes the quote when passed to powershell.exe. Strip it so paths with
rem spaces are quoted correctly.
if "%PROJECT_DIR:~-1%"=="\" set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"
set "LAUNCHER=%PROJECT_DIR%\scripts\run_cellpose_gui_windows.ps1"

if not exist "%LAUNCHER%" (
  echo The Windows launcher helper was not found:
  echo %LAUNCHER%
  echo.
  pause
  exit /b 1
)

powershell.exe -NoProfile -STA -ExecutionPolicy Bypass -File "%LAUNCHER%" -ProjectDir "%PROJECT_DIR%"
set "STATUS=%ERRORLEVEL%"

if not "%STATUS%"=="0" (
  echo.
  echo The Cellpose DAPI / aSMA app did not start successfully.
  echo.
  pause
)

exit /b %STATUS%
