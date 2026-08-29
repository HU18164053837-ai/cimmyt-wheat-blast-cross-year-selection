$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$required = @('README.md','LICENSE','DATA_LICENSE.md','CITATION.cff','requirements.txt','RUN_ALL.ps1','RUN_ORDER.md','PARAMETER_DICTIONARY.md')
$failed = $false

foreach ($name in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $root $name))) {
        Write-Error "Missing required file: $name"
        $failed = $true
    }
}

$files = Get-ChildItem -LiteralPath $root -Recurse -File | Where-Object { $_.FullName -notmatch '[\\/]\.git[\\/]' }
$oversized = $files | Where-Object Length -ge 95MB
if ($oversized) {
    $oversized | ForEach-Object { Write-Error "File approaches/exceeds GitHub limit: $($_.FullName)" }
    $failed = $true
}

$textFiles = $files | Where-Object { $_.Name -ne 'CHECK_BEFORE_UPLOAD.ps1' -and $_.Extension -in @('.R','.py','.ps1','.md','.txt','.csv','.yml','.yaml','.cff') -and $_.Length -lt 20MB }
$patterns = @(
    'C:\\Users\\',
    '/Users/[^/]+/',
    '18164053837',
    '(?i)(password|passwd|api[_-]?key|secret)[ ]*[:=][ ]*["''][^"'']+["'']',
    'BEGIN (RSA|OPENSSH|EC|DSA) PRIVATE KEY'
)
foreach ($pattern in $patterns) {
    $hits = $textFiles | Select-String -Pattern $pattern -ErrorAction SilentlyContinue
    if ($hits) {
        $hits | ForEach-Object { Write-Error "Potential private/local content: $($_.Path):$($_.LineNumber)" }
        $failed = $true
    }
}

$inventory = foreach ($file in $files) {
    [pscustomobject]@{
        path = $file.FullName.Substring($root.Length + 1).Replace('\\','/')
        bytes = $file.Length
        sha256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    }
}
$inventory | Export-Csv -LiteralPath (Join-Path $root 'REPOSITORY_FILE_MANIFEST.csv') -NoTypeInformation -Encoding UTF8

if ($failed) { exit 1 }
Write-Host ("Pre-upload checks passed: {0} files, {1:N2} MB" -f $files.Count,(($files | Measure-Object Length -Sum).Sum / 1MB)) -ForegroundColor Green
exit 0
