param(
    [bool]$InstallPythonPackages = $true
)

$ErrorActionPreference = 'Stop'

function Write-Section([string]$title) {
    Write-Host "`n=== $title ===" -ForegroundColor Cyan
}

function Test-CommandAvailable([string]$name, [switch]$Required) {
    $cmd = Get-Command $name -ErrorAction SilentlyContinue
    if ($null -eq $cmd) {
        if ($Required) {
            throw "Required command '$name' is missing from PATH."
        }
        Write-Host "MISSING: $name" -ForegroundColor Yellow
        return $null
    }
    Write-Host "FOUND: $name -> $($cmd.Source)" -ForegroundColor Green
    return $cmd.Source
}

Write-Section "Checking command availability on PATH"
Test-CommandAvailable -name 'python' -Required | Out-Null
Test-CommandAvailable -name 'java' -Required | Out-Null
Test-CommandAvailable -name 'quarto' -Required | Out-Null

$windowsAppsPath = Join-Path $env:LocalAppData 'Microsoft\WindowsApps'

# These are required for this repo workflows but may not be present initially
$jupyterPath = Test-CommandAvailable -name 'jupyter'
$mvnPath = Test-CommandAvailable -name 'mvn'

if ($InstallPythonPackages) {
    Write-Section "Installing/upgrading notebook Python packages"
    $reqFile = Join-Path $PSScriptRoot 'requirements-fiba-notebook.txt'
    if (-not (Test-Path $reqFile)) {
        throw "Missing requirements file: $reqFile"
    }

    & python -m pip install --upgrade pip
    & python -m pip install --upgrade -r $reqFile

    Write-Section "Ensuring user Scripts directory is on PATH"
    $scriptsDir = & python -c "import site; print(site.USER_BASE + r'\\Scripts')"
    if (-not (Test-Path $scriptsDir)) {
        New-Item -ItemType Directory -Path $scriptsDir | Out-Null
    }

    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if ([string]::IsNullOrWhiteSpace($userPath)) { $userPath = '' }

    $normalized = ($userPath -split ';' | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' })
    if ($normalized -notcontains $scriptsDir) {
        $newPath = if ($userPath.Trim().Length -gt 0) { "$userPath;$scriptsDir" } else { $scriptsDir }
        [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')
        Write-Host "Added to USER PATH: $scriptsDir" -ForegroundColor Green
        Write-Host "Restart VS Code/terminal to load updated PATH." -ForegroundColor Yellow
    }
    else {
        Write-Host "USER PATH already contains: $scriptsDir" -ForegroundColor Green
    }

    Write-Section "Removing WindowsApps shims from USER PATH (default)"
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if ([string]::IsNullOrWhiteSpace($userPath)) { $userPath = '' }

    $parts = @(
        $userPath -split ';' |
            ForEach-Object { $_.Trim() } |
            Where-Object { $_ -and ($_ -ne $windowsAppsPath) }
    )

    $newUserPath = ($parts -join ';')
    if ($newUserPath -ne $userPath) {
        [Environment]::SetEnvironmentVariable('Path', $newUserPath, 'User')
        Write-Host "Removed from USER PATH: $windowsAppsPath" -ForegroundColor Green
        Write-Host "Restart VS Code/terminal to load updated PATH." -ForegroundColor Yellow
    }
    else {
        Write-Host "USER PATH already excludes: $windowsAppsPath" -ForegroundColor Green
    }
}

Write-Section "Post-check"
$jupyterPath = Test-CommandAvailable -name 'jupyter'
$mvnPath = Test-CommandAvailable -name 'mvn'

if ($null -eq $jupyterPath) {
    Write-Host "jupyter still not found on PATH; restart shell, then rerun this script." -ForegroundColor Yellow
}
if ($null -eq $mvnPath) {
    Write-Host "mvn is still missing. Install Maven and add its bin/ directory to PATH." -ForegroundColor Yellow
}

Write-Host "`nEnvironment setup script completed." -ForegroundColor Cyan
