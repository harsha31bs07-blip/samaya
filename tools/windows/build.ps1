# Builds samaya with MSVC: finds Visual Studio (vswhere), enters its developer shell (cl, CMake,
# Ninja on PATH), then configures, builds and optionally tests one CMake preset.
#
#   powershell -ExecutionPolicy Bypass -File tools\windows\build.ps1 [-Preset windows-release] [-Test]
param(
  [string]$Preset = 'windows-release',
  [switch]$Test,
  [switch]$KeepGoing,  # build every target despite errors (to see all diagnostics at once)
  [string]$Source = (Split-Path -Parent (Split-Path -Parent $PSScriptRoot))
)
$ErrorActionPreference = 'Stop'

$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path $vswhere)) { throw 'Visual Studio Build Tools are not installed (run setup.bat).' }
$vs = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
  -property installationPath
if (-not $vs) { throw 'The MSVC C++ tools are not installed (run setup.bat).' }

Import-Module (Join-Path $vs 'Common7\Tools\Microsoft.VisualStudio.DevShell.dll')
Enter-VsDevShell -VsInstallPath $vs -SkipAutomaticLocation -DevCmdArguments '-arch=x64 -host_arch=x64' |
  Out-Null

Set-Location $Source
cmake --preset $Preset
if ($LASTEXITCODE) { exit $LASTEXITCODE }
if ($KeepGoing) { cmake --build --preset $Preset -- -k 0 } else { cmake --build --preset $Preset }
if ($LASTEXITCODE) { exit $LASTEXITCODE }
if ($Test) {
  ctest --preset $Preset
  exit $LASTEXITCODE
}
exit 0
