@echo off
setlocal

set "PROJECT_DIR=%~dp0"
rem %~dp0 always ends in a backslash; a trailing "\" before the closing quote
rem escapes the quote when passed to powershell.exe. Strip it so paths with
rem spaces are quoted correctly.
if "%PROJECT_DIR:~-1%"=="\" set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"
set "SETUP_SCRIPT=%PROJECT_DIR%\scripts\setup_cellpose_pipeline_windows.ps1"

if not exist "%SETUP_SCRIPT%" (
  echo The Windows setup helper was not found:
  echo %SETUP_SCRIPT%
  echo.
  pause
  exit /b 1
)

echo Project: %PROJECT_DIR%
echo.
echo Setting up the project-local Cellpose environment. This can take a while.
echo.

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%SETUP_SCRIPT%" -ProjectDir "%PROJECT_DIR%"
set "STATUS=%ERRORLEVEL%"

if not "%STATUS%"=="0" (
  echo.
  echo Setup failed. The Command Prompt window contains the error details.
  echo.
  pause
)

exit /b %STATUS%
