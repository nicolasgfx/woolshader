@echo off
REM -------------------------------------------------------------------
REM  build.bat -- configure (if needed) and build Release with MSVC + Ninja
REM
REM    build.bat                    build everything
REM    build.bat test               build + run ctest
REM    build.bat check              build + run knit_check --all
REM    build.bat run                build + launch the viewer
REM    build.bat nogui              configure without the viewer (no downloads)
REM    build.bat build <args>       build with extra cmake --build arguments, e.g.
REM                                 build.bat build --target knit_tests_ui
REM
REM  Environment:
REM    KNIT_BUILD_DIR   build folder (default: <checkout>-build\release next to the checkout). A folder configured from another
REM                     checkout is refused, so worktrees never share a build.
REM    KNIT_CUDA_ROOT   CUDA toolkit (default: CUDA v13.4)
REM    KNIT_CMAKE_ARGS  extra cmake options when a build folder is configured
REM
REM  The build includes the CUDA backend (GPU knitting, GPU fibres).
REM -------------------------------------------------------------------
setlocal enabledelayedexpansion
cd /d "%~dp0"

set "BUILD_DIR=%KNIT_BUILD_DIR%"
if "%BUILD_DIR%"=="" for %%I in ("%~dp0.") do set "BUILD_DIR=%%~dpI%%~nxI-build\release"
set "MODE=%~1"
set CMAKE_EXTRA=
if /I "%MODE%"=="nogui" set CMAKE_EXTRA=-DKNIT_BUILD_VIEWER=OFF
REM  KNIT_CMAKE_ARGS: extra configure options for a new build folder, e.g. -DKNIT_BUILD_TESTS=OFF
if defined KNIT_CMAKE_ARGS set CMAKE_EXTRA=%CMAKE_EXTRA% %KNIT_CMAKE_ARGS%
set "SRC=%~dp0"
set BUILD_ARGS=
if /I not "%MODE%"=="build" goto args_done
shift /1
:collect
if "%~1"=="" goto args_done
set BUILD_ARGS=!BUILD_ARGS! %1
shift /1
goto collect
:args_done

REM -- the build folder belongs to this checkout -----------------------
set "SRC=%SRC:~0,-1%"
set "SRC=%SRC:\=/%"
if exist "%BUILD_DIR%\CMakeCache.txt" (
    for /f "tokens=2 delims==" %%D in ('findstr /B /C:"CMAKE_HOME_DIRECTORY:INTERNAL=" "%BUILD_DIR%\CMakeCache.txt"') do set "HOME_DIR=%%D"
    if /I not "!HOME_DIR!"=="%SRC%" (
        echo [build] %BUILD_DIR% was configured for !HOME_DIR!, not this checkout ^(%SRC%^): set KNIT_BUILD_DIR
        exit /b 1
    )
)

REM -- CUDA toolkit -----------------------------------------------------
set "CUDA=%KNIT_CUDA_ROOT%"
if "%CUDA%"=="" set "CUDA=%ProgramFiles%\NVIDIA GPU Computing Toolkit\CUDA\v13.4"
if not exist "%CUDA%\bin\nvcc.exe" (
    echo [build] CUDA toolkit not found at "%CUDA%" ^(set KNIT_CUDA_ROOT^)
    exit /b 1
)
set CMAKE_GPU="-DCMAKE_CUDA_COMPILER=%CUDA%\bin\nvcc.exe" "-DCUDAToolkit_ROOT=%CUDA%"

REM -- MSVC environment ------------------------------------------------
REM  vcvarsall.bat itself calls vswhere.exe without a path; put the VS Installer
REM  directory on PATH first, otherwise it prints "vswhere.exe not found".
set "PATH=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer;%PATH%"
where cl.exe >nul 2>&1
if errorlevel 1 (
    set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
    if not exist "!VSWHERE!" ( echo [build] vswhere not found & exit /b 1 )
    for /f "usebackq tokens=*" %%I in (`"!VSWHERE!" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSROOT=%%I"
    if not defined VSROOT ( echo [build] No VS C++ toolset found & exit /b 1 )
    call "!VSROOT!\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul
)
REM -- VS-bundled CMake / Ninja -----------------------------------------
where cmake.exe >nul 2>&1
if errorlevel 1 (
    if not defined VSROOT for /f "usebackq tokens=*" %%I in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -property installationPath`) do set "VSROOT=%%I"
    set "PATH=!VSROOT!\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin;!VSROOT!\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja;!PATH!"
)

set RECONFIGURE=
if /I "%MODE%"=="nogui" set RECONFIGURE=1
if not exist "%BUILD_DIR%\build.ninja" (
    cmake -S . -B "%BUILD_DIR%" -G Ninja -DCMAKE_BUILD_TYPE=Release %CMAKE_GPU% %CMAKE_EXTRA%
    if errorlevel 1 exit /b 1
) else if defined RECONFIGURE (
    echo [build] configuring %BUILD_DIR%
    cmake -S . -B "%BUILD_DIR%" %CMAKE_GPU% %CMAKE_EXTRA%
    if errorlevel 1 exit /b 1
)

cmake --build "%BUILD_DIR%" -j %NUMBER_OF_PROCESSORS% %BUILD_ARGS%
if errorlevel 1 exit /b 1

if /I "%MODE%"=="test" (
    ctest --test-dir "%BUILD_DIR%" --output-on-failure
    if errorlevel 1 exit /b 1
)
if /I "%MODE%"=="check" (
    "%BUILD_DIR%\knit_check.exe" --all
    if errorlevel 1 exit /b 1
)
if /I "%MODE%"=="run" (
    start "" "%BUILD_DIR%\knit_viewer.exe"
)
echo [build] OK
endlocal
