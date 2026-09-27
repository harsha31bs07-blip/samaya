@echo off
rem samaya setup: checks this PC for what samaya needs, installs what is missing (asking first),
rem then installs samaya Studio. Double-click it, or run:  setup.bat [/check] [/yes] [/uninstall]
setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\windows\setup.ps1" %*
set RESULT=%ERRORLEVEL%
echo.
if not "%RESULT%"=="0" echo Setup did not finish; the messages above say why. The log is in %LOCALAPPDATA%\samaya\setup.log.
pause
exit /b %RESULT%
