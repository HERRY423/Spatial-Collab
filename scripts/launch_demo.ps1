$ErrorActionPreference = 'Stop'
$pluginRoot = Split-Path -Parent $PSScriptRoot
$env:PYTHONPATH = Join-Path $pluginRoot 'src'
$env:PYTHONIOENCODING = 'utf-8'
$projectPath = Join-Path $pluginRoot 'projects\demo'
if (-not (Test-Path -LiteralPath $projectPath)) {
    python -m spatial_collab demo $projectPath
    if ($LASTEXITCODE -ne 0) { throw 'Demo creation failed.' }
}
Write-Host 'Open http://127.0.0.1:8765 in your browser. Keep this process running while reviewing.'
python -m spatial_collab serve --project $projectPath --port 8765

