@echo off
REM -------------------------------------------------------------------
REM  run.bat -- build (Release) and launch the knit viewer
REM
REM    run.bat                      build + start the viewer
REM    run.bat --scene teddy        any arguments are passed to knit_viewer
REM  KNIT_BUILD_DIR as in build.bat (default: <checkout>-build\release next to the checkout).
REM -------------------------------------------------------------------
setlocal
cd /d "%~dp0"
set "BUILD_DIR=%KNIT_BUILD_DIR%"
if "%BUILD_DIR%"=="" for %%I in ("%~dp0.") do set "BUILD_DIR=%%~dpI%%~nxI-build\release"

call "%~dp0build.bat"
if errorlevel 1 (
    echo [run] build failed
    exit /b 1
)

if not exist "%BUILD_DIR%\knit_viewer.exe" (
    echo [run] %BUILD_DIR%\knit_viewer.exe not found -- was the build configured with "nogui"?
    echo [run] delete the build folder and run again to re-enable the viewer
    exit /b 1
)

"%BUILD_DIR%\knit_viewer.exe" %*
endlocal
