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

# The CUDA Toolkit (for the windows-cuda preset). A shell started before it was installed lacks its
# variables, so take them from the machine environment.
if (-not $env:CUDA_PATH) { $env:CUDA_PATH = [Environment]::GetEnvironmentVariable('CUDA_PATH', 'Machine') }
if ($env:CUDA_PATH) { $env:PATH = (Join-Path $env:CUDA_PATH 'bin') + ';' + $env:PATH }
# nvcc fails without a message when TEMP has a space in it (a user folder such as C:\Users\A B):
# use the folder's short name, or else a folder at the root of the system drive.
if ($env:TEMP -match ' ') {
  $short = (New-Object -ComObject Scripting.FileSystemObject).GetFolder($env:TEMP).ShortPath
  if ($short -match ' ') { $short = Join-Path $env:SystemDrive 'samaya-tmp'; New-Item -ItemType Directory -Force $short | Out-Null }
  $env:TEMP = $short
  $env:TMP = $short
}

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
