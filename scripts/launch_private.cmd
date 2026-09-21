@echo off
setlocal
if not defined SPATIAL_COLLAB_HOME set "SPATIAL_COLLAB_HOME=%LOCALAPPDATA%\SpatialCollab"
set "SPATIAL_PRIVATE_PYTHON=%SPATIAL_COLLAB_HOME%\runtime\Scripts\python.exe"
if not exist "%SPATIAL_PRIVATE_PYTHON%" (
  >&2 echo Spatial Collab is not configured. Run python scripts/setup_private.py --create-empty --with-import first.
  exit /b 2
)
set "PYTHONPATH=%~dp0..\src"
set "PYTHONIOENCODING=utf-8"
set "SPATIAL_PRIVATE_ACTION=mcp"
if "%~1"=="--doctor" set "SPATIAL_PRIVATE_ACTION=doctor"
"%SPATIAL_PRIVATE_PYTHON%" -m spatial_collab.private_runtime %SPATIAL_PRIVATE_ACTION%
exit /b %errorlevel%
