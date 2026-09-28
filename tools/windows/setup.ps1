# samaya setup for Windows (run through setup.bat). Checks this PC for what samaya needs,
# installs what is missing (asking first), then installs samaya Studio for the current user.
#
# Two modes, chosen by what sits next to setup.bat:
#   release package (bin\samaya-studio.exe): installs the prebuilt app; needs only the WebView2
#     Runtime (part of Windows 11 and current Windows 10).
#   source tree (CMakeLists.txt): also needs Visual Studio 2022 Build Tools (C++), and Python 3
#     for the sample cases; builds with MSVC, runs the tests, then installs.
#
#   setup.bat [/check] [/yes] [/uninstall] [/nolaunch]
#     /check      report what is present and missing; change nothing
#     /yes        install missing components without asking
#     /uninstall  remove samaya Studio (this user)
#     /nolaunch   do not open Studio at the end
# Windows PowerShell 5.1 compatible. Logs to %LOCALAPPDATA%\samaya\setup.log.

$ErrorActionPreference = 'Stop'
$opts = @{}
foreach ($a in $args) { $opts[($a -replace '^[/-]+', '').ToLower()] = $true }
$CheckOnly = $opts.ContainsKey('check')
$AssumeYes = $opts.ContainsKey('yes')

$Root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$InstallDir = Join-Path $env:LOCALAPPDATA 'Programs\samaya'
$DataDir = Join-Path $env:LOCALAPPDATA 'samaya'
New-Item -ItemType Directory -Force -Path $DataDir | Out-Null
try { Start-Transcript -Path (Join-Path $DataDir 'setup.log') -Append | Out-Null } catch { }

function Say($text) { Write-Host $text }
function Ok($text) { Write-Host ('  [ok]      ' + $text) -ForegroundColor Green }
function Missing($text) { Write-Host ('  [missing] ' + $text) -ForegroundColor Yellow }
function Bad($text) { Write-Host ('  [error]   ' + $text) -ForegroundColor Red }
function Ask($question) {
  if ($AssumeYes) { return $true }
  $answer = Read-Host ($question + ' [Y/n]')
  return ($answer -eq '' -or $answer -match '^[yY]')
}

# --- Checks ---------------------------------------------------------------------------------------
function Get-WebView2Version {
  $keys = @('HKLM:\SOFTWARE\WOW6432Node\Microsoft\EdgeUpdate\Clients\{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}',
            'HKLM:\SOFTWARE\Microsoft\EdgeUpdate\Clients\{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}',
            'HKCU:\Software\Microsoft\EdgeUpdate\Clients\{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}')
  foreach ($k in $keys) {
    $pv = (Get-ItemProperty -Path $k -Name pv -ErrorAction SilentlyContinue).pv
    if ($pv -and $pv -ne '0.0.0.0') { return $pv }
  }
  return $null
}
function Get-VsPath {
  $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
  if (-not (Test-Path $vswhere)) { return $null }
  $p = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
  if ($p) { return $p } else { return $null }
}
function Get-Python {
  foreach ($cmd in @('py', 'python', 'python3')) {
    $c = Get-Command $cmd -ErrorAction SilentlyContinue
    if (-not $c) { continue }
    try {
      $args3 = @(); if ($cmd -eq 'py') { $args3 = @('-3') }
      $v = & $cmd @args3 --version 2>&1
      if ($LASTEXITCODE -eq 0 -and "$v" -match 'Python 3\.(\d+)' -and [int]$Matches[1] -ge 9) {
        return @{ Cmd = $cmd; Args = $args3; Version = "$v".Trim() }
      }
    } catch { }
  }
  return $null
}
function Get-NvidiaGpu {
  try {
    $g = Get-CimInstance Win32_VideoController | Where-Object { $_.Name -match 'NVIDIA' } | Select-Object -First 1
    if ($g) { return $g.Name }
  } catch { }
  return $null
}
function Refresh-Path {
  $env:Path = [Environment]::GetEnvironmentVariable('Path', 'Machine') + ';' + [Environment]::GetEnvironmentVariable('Path', 'User')
}
function Install-WithWinget($id, $override) {
  if (-not (Get-Command winget -ErrorAction SilentlyContinue)) { throw 'winget (App Installer) is not available on this PC.' }
  $a = @('install', '--id', $id, '-e', '--accept-package-agreements', '--accept-source-agreements')
  if ($override) { $a += @('--override', $override) }
  & winget @a
  # 3010: installed, a restart is needed later.
  if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne 3010 -and $LASTEXITCODE -ne -1978335189) {
    throw "winget could not install $id (exit code $LASTEXITCODE)."
  }
  Refresh-Path
}
function Install-WebView2 {
  if (Get-Command winget -ErrorAction SilentlyContinue) {
    try { Install-WithWinget 'Microsoft.EdgeWebView2Runtime' $null; return } catch { Say ('  winget: ' + $_.Exception.Message) }
  }
  # Microsoft's Evergreen bootstrapper, checked to be signed by Microsoft before it runs.
  $exe = Join-Path $env:TEMP 'MicrosoftEdgeWebview2Setup.exe'
  Invoke-WebRequest -UseBasicParsing -Uri 'https://go.microsoft.com/fwlink/p/?LinkId=2124703' -OutFile $exe
  $sig = Get-AuthenticodeSignature $exe
  if ($sig.Status -ne 'Valid' -or $sig.SignerCertificate.Subject -notmatch 'Microsoft') { throw 'The WebView2 installer is not signed by Microsoft; not running it.' }
  $p = Start-Process -FilePath $exe -ArgumentList '/silent', '/install' -Wait -PassThru
  if ($p.ExitCode -ne 0) { throw "The WebView2 installer failed (exit code $($p.ExitCode))." }
}

# --- Install and uninstall ------------------------------------------------------------------------
function New-Shortcut($path, $target, $workdir) {
  $shell = New-Object -ComObject WScript.Shell
  $s = $shell.CreateShortcut($path)
  $s.TargetPath = $target
  $s.WorkingDirectory = $workdir
  $s.IconLocation = "$target,0"
  $s.Description = 'samaya Studio: verified LP, MILP and QP solver'
  $s.Save()
}
function Install-Studio($from) {
  Say ''
  Say "Installing samaya Studio to $InstallDir"
  Get-Process samaya-studio -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
  New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
  foreach ($item in @('samaya.exe', 'samaya-studio.exe')) { Copy-Item -Force (Join-Path $from $item) $InstallDir }
  foreach ($dir in @('ui', 'samples')) {
    $src = Join-Path $from $dir
    if (Test-Path $src) {
      Remove-Item -Recurse -Force (Join-Path $InstallDir $dir) -ErrorAction SilentlyContinue
      Copy-Item -Recurse -Force $src (Join-Path $InstallDir $dir)
    }
  }
  $exe = Join-Path $InstallDir 'samaya-studio.exe'
  $startMenu = Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs'
  New-Shortcut (Join-Path $startMenu 'samaya Studio.lnk') $exe $InstallDir
  New-Shortcut (Join-Path ([Environment]::GetFolderPath('Desktop')) 'samaya Studio.lnk') $exe $InstallDir
  Ok 'samaya Studio installed (Start menu and desktop shortcuts)'
  Ok ("command-line solver: " + (Join-Path $InstallDir 'samaya.exe'))
}
function Uninstall-Studio {
  Get-Process samaya-studio -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
  Remove-Item -Recurse -Force $InstallDir -ErrorAction SilentlyContinue
  Remove-Item -Force (Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs\samaya Studio.lnk') -ErrorAction SilentlyContinue
  Remove-Item -Force (Join-Path ([Environment]::GetFolderPath('Desktop')) 'samaya Studio.lnk') -ErrorAction SilentlyContinue
  Ok 'samaya Studio removed (your run history in %LOCALAPPDATA%\samaya is kept)'
}

# --- Main -------------------------------------------------------------------------------------------
try {
  Say ''
  Say 'samaya setup'
  Say '============'
  if ($opts.ContainsKey('uninstall')) { Uninstall-Studio; exit 0 }

  $release = Test-Path (Join-Path $Root 'bin\samaya-studio.exe')
  $source = (-not $release) -and (Test-Path (Join-Path $Root 'CMakeLists.txt'))
  if (-not $release -and -not $source) { throw "Run setup.bat from the samaya release folder or the source tree ($Root)." }
  if ($release) { Say 'Mode: release package (prebuilt)' } else { Say 'Mode: source tree (build with MSVC)' }
  Say ''
  Say 'Checking this PC'

  # Windows version and architecture.
  $os = [Environment]::OSVersion.Version
  if (-not [Environment]::Is64BitOperatingSystem) { Bad '64-bit Windows is required'; exit 1 }
  if ($os.Major -lt 10 -or ($os.Major -eq 10 -and $os.Build -lt 17763)) { Bad "Windows 10 (1809) or later is required; this is $os"; exit 1 }
  $winName = if ($os.Build -ge 22000) { 'Windows 11' } else { 'Windows 10' }  # 11 reports 10.0
  Ok "$winName (build $($os.Build)), 64-bit"

  $todo = @()
  $wv = Get-WebView2Version
  if ($wv) { Ok "Microsoft Edge WebView2 Runtime $wv" } else { Missing 'Microsoft Edge WebView2 Runtime (needed by samaya Studio; about 2 MB download)'; $todo += 'webview2' }

  $py = $null
  if ($source) {
    if (Get-Command winget -ErrorAction SilentlyContinue) { Ok 'winget (App Installer)' } else { Missing 'winget (App Installer): install "App Installer" from the Microsoft Store' }
    $vs = Get-VsPath
    if ($vs) { Ok "Visual Studio C++ Build Tools ($vs)" } else { Missing 'Visual Studio 2022 Build Tools with C++ (compiler, Windows SDK, CMake, Ninja; several GB)'; $todo += 'vs' }
    $py = Get-Python
    if ($py) { Ok "$($py.Version) (sample cases and reports)" } else { Missing 'Python 3.9 or later (sample cases and reports; about 100 MB)'; $todo += 'python' }
  }
  $gpu = Get-NvidiaGpu
  if ($gpu) { Ok "GPU: $gpu (samaya Studio can run PDLP on it: Settings, Compute)" } else { Say '  [info]    no NVIDIA GPU found; samaya runs on the CPU' }

  if ($CheckOnly) {
    Say ''
    if ($todo.Count -eq 0) { Say 'Everything samaya needs is present.' } else { Say ("Missing: " + ($todo -join ', ') + ". Run setup.bat without /check to install.") }
    exit 0
  }

  foreach ($t in $todo) {
    Say ''
    if ($t -eq 'webview2' -and (Ask 'Install the Microsoft Edge WebView2 Runtime now?')) { Install-WebView2; Ok 'WebView2 Runtime installed' }
    elseif ($t -eq 'vs' -and (Ask 'Install Visual Studio 2022 Build Tools (C++) now? This downloads several GB and asks for administrator approval.')) {
      Install-WithWinget 'Microsoft.VisualStudio.2022.BuildTools' '--quiet --wait --norestart --nocache --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended --add Microsoft.VisualStudio.Component.VC.CMake.Project --add Microsoft.VisualStudio.Component.VC.ASAN'
      if (-not (Get-VsPath)) { throw 'The Build Tools installation did not finish; run setup.bat again (a restart may be needed).' }
      Ok 'Visual Studio Build Tools installed'
    }
    elseif ($t -eq 'python' -and (Ask 'Install Python 3.12 now?')) { Install-WithWinget 'Python.Python.3.12' $null; $py = Get-Python; if ($py) { Ok $py.Version } }
    elseif ($t -eq 'webview2' -or $t -eq 'vs') { throw "$t is required; setup stopped." }
  }

  if ($release) {
    Install-Studio (Join-Path $Root 'bin')
  } else {
    Say ''
    Say 'Building samaya with MSVC (release) and running the tests'
    & powershell -NoProfile -ExecutionPolicy Bypass -File (Join-Path $PSScriptRoot 'build.ps1') -Preset windows-release -Test
    if ($LASTEXITCODE -ne 0) { throw 'The build or its tests failed; see above.' }
    Ok 'build and tests passed'
    $out = Join-Path $Root 'build\windows-release\apps\studio'
    if ($py) {
      $samples = Join-Path $out 'samples'
      New-Item -ItemType Directory -Force -Path $samples | Out-Null
      & $py.Cmd @($py.Args) (Join-Path $Root 'cases\mrpl.py') generate --out $samples | Out-Null
      Ok 'sample refinery cases generated'
    }
    Install-Studio $out
  }

  if (-not $opts.ContainsKey('nolaunch')) { Start-Process (Join-Path $InstallDir 'samaya-studio.exe') }
  Say ''
  Say 'Done. Open "samaya Studio" from the Start menu or the desktop.'
  exit 0
} catch {
  Bad $_.Exception.Message
  exit 1
} finally {
  try { Stop-Transcript | Out-Null } catch { }
}
