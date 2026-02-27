param(
    [string]$PythonRoot = '',
    [string]$MavenBin = '',
    [switch]$WhatIf
)

$ErrorActionPreference = 'Stop'

function Write-Section([string]$Title) {
    Write-Host "`n=== $Title ===" -ForegroundColor Cyan
}

function Get-CleanPathEntry([string]$Entry) {
    if ([string]::IsNullOrWhiteSpace($Entry)) { return $null }
    $v = $Entry.Trim().Trim('"')
    if ($v.Length -eq 0) { return $null }
    return $v
}

function Add-UniquePathEntry {
    param(
        [System.Collections.Generic.List[string]]$List,
        [string]$Entry,
        [bool]$MustExist = $false,
        [string]$Reason = ''
    )

    $normalized = Get-CleanPathEntry $Entry
    if ($null -eq $normalized) { return }

    if ($MustExist -and -not (Test-Path -LiteralPath $normalized)) {
        Write-Host "SKIP (missing): $normalized" -ForegroundColor Yellow
        return
    }

    $already = $false
    foreach ($existing in $List) {
        if ($existing.Equals($normalized, [System.StringComparison]::OrdinalIgnoreCase)) {
            $already = $true
            break
        }
    }

    if (-not $already) {
        $List.Add($normalized)
        if ($Reason) {
            Write-Host "ADD   : $normalized  [$Reason]" -ForegroundColor Green
        }
        else {
            Write-Host "ADD   : $normalized" -ForegroundColor Green
        }
    }
    else {
        if ($Reason) {
            Write-Host "KEEP  : $normalized  [$Reason]" -ForegroundColor DarkGray
        }
    }
}

function Find-PythonRootAuto {
    $candidates = @(
        "$env:LocalAppData\\Programs\\Python",
        "$env:ProgramFiles",
        "${env:ProgramFiles(x86)}"
    )

    $best = $null
    foreach ($base in $candidates) {
        if (-not (Test-Path -LiteralPath $base)) { continue }

        $dirs = Get-ChildItem -LiteralPath $base -Directory -ErrorAction SilentlyContinue |
            Where-Object { $_.Name -match '^Python\d{2,3}$|^Python\s*\d' }

        foreach ($d in $dirs) {
            $py = Join-Path $d.FullName 'python.exe'
            if (Test-Path -LiteralPath $py) {
                if ($null -eq $best) { $best = $d.FullName }
                else {
                    # Prefer lexical max folder name (good enough for Python311 vs Python312)
                    if ($d.FullName -gt $best) { $best = $d.FullName }
                }
            }
        }
    }
    return $best
}

function Find-MavenBinAuto {
    $bases = @(
        "$env:ProgramFiles",
        "${env:ProgramFiles(x86)}",
        "$env:LocalAppData\\Programs",
        "C:\\Tools",
        "C:\\dev"
    )

    foreach ($base in $bases) {
        if (-not (Test-Path -LiteralPath $base)) { continue }

        $dirs = Get-ChildItem -LiteralPath $base -Directory -ErrorAction SilentlyContinue |
            Where-Object { $_.Name -match 'apache-maven|maven' }

        foreach ($d in $dirs) {
            $bin = Join-Path $d.FullName 'bin'
            $mvnCmd = Join-Path $bin 'mvn.cmd'
            $mvnExe = Join-Path $bin 'mvn.exe'
            if ((Test-Path -LiteralPath $mvnCmd) -or (Test-Path -LiteralPath $mvnExe)) {
                return $bin
            }
        }
    }
    return $null
}

Write-Section 'Read current USER PATH'
$currentUserPath = [Environment]::GetEnvironmentVariable('Path', 'User')
if ($null -eq $currentUserPath) { $currentUserPath = '' }
Write-Host "Current USER PATH length: $($currentUserPath.Length)"

$currentEntries = New-Object 'System.Collections.Generic.List[string]'
foreach ($raw in ($currentUserPath -split ';')) {
    $n = Get-CleanPathEntry $raw
    if ($null -ne $n) { $currentEntries.Add($n) }
}

# Resolve Python root
if ([string]::IsNullOrWhiteSpace($PythonRoot)) {
    $PythonRoot = Find-PythonRootAuto
}
if (-not [string]::IsNullOrWhiteSpace($PythonRoot)) {
    $PythonRoot = (Get-CleanPathEntry $PythonRoot)
}
$PythonScripts = if ($PythonRoot) { Join-Path $PythonRoot 'Scripts' } else { $null }

# Resolve Maven bin
if ([string]::IsNullOrWhiteSpace($MavenBin)) {
    $MavenBin = Find-MavenBinAuto
}
if (-not [string]::IsNullOrWhiteSpace($MavenBin)) {
    $MavenBin = (Get-CleanPathEntry $MavenBin)
}

Write-Section 'Build surgically reordered PATH'
$newEntries = New-Object 'System.Collections.Generic.List[string]'

# 1) Put real Python first (if installed)
if ($PythonRoot) {
    Add-UniquePathEntry -List $newEntries -Entry $PythonRoot -MustExist $true -Reason 'real python root first'
    Add-UniquePathEntry -List $newEntries -Entry $PythonScripts -MustExist $true -Reason 'python Scripts for pip/jupyter'
}
else {
    Write-Host 'No real Python installation detected automatically. Provide -PythonRoot explicitly after install.' -ForegroundColor Yellow
}

# 2) Keep good known tools near top
Add-UniquePathEntry -List $newEntries -Entry "$env:LocalAppData\\Programs\\Eclipse Adoptium\\jdk-25.0.2.10-hotspot\\bin" -MustExist $false -Reason 'java'
Add-UniquePathEntry -List $newEntries -Entry "$env:LocalAppData\\Programs\\Quarto\\bin" -MustExist $false -Reason 'quarto'

# 3) Add Maven bin if available
if ($MavenBin) {
    Add-UniquePathEntry -List $newEntries -Entry $MavenBin -MustExist $true -Reason 'maven'
}
else {
    Write-Host 'No Maven bin detected automatically. Provide -MavenBin explicitly after install.' -ForegroundColor Yellow
}

# 4) Preserve everything else except WindowsApps
$windowsApps = "$env:LocalAppData\\Microsoft\\WindowsApps"
foreach ($entry in $currentEntries) {
    if ($entry.Equals($windowsApps, [System.StringComparison]::OrdinalIgnoreCase)) {
        continue
    }
    Add-UniquePathEntry -List $newEntries -Entry $entry -MustExist $false
}

# 5) Deliberately do NOT add WindowsApps back

$newUserPath = [string]::Join(';', $newEntries)

Write-Section 'Diff preview'
Write-Host '--- OLD USER PATH ---'
Write-Host $currentUserPath
Write-Host "`n--- NEW USER PATH ---"
Write-Host $newUserPath

if ($WhatIf) {
    Write-Section 'WhatIf mode'
    Write-Host 'No changes written. Re-run without -WhatIf to apply.' -ForegroundColor Yellow
    exit 0
}

Write-Section 'Apply USER PATH update'
[Environment]::SetEnvironmentVariable('Path', $newUserPath, 'User')
Write-Host 'USER PATH updated successfully.' -ForegroundColor Green

Write-Section 'Next'
Write-Host '1) Restart VS Code and terminals.'
Write-Host '2) Verify: python, pip, jupyter, mvn, java, quarto'
Write-Host '3) If python still resolves to WindowsApps, disable Python App Execution Aliases in Windows Settings.'
Write-Host '4) This script intentionally excludes WindowsApps from USER PATH.'
