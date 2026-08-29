param(
    [string]$RepositoryName = 'cimmyt-wheat-blast-cross-year-selection',
    [ValidateSet('public','private')][string]$Visibility = 'public',
    [string]$Description = 'Code and reproducibility materials for auditable cross-year selection of wheat blast resistance candidates'
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $root

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Git is not installed or is not available on PATH.'
}
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw 'GitHub CLI (gh) is not installed. Install it from https://cli.github.com/ and run gh auth login.'
}

gh auth status
if ($LASTEXITCODE -ne 0) { throw 'GitHub CLI is not logged in. Run: gh auth login' }
$owner = (gh api user --jq '.login').Trim()
if (-not $owner) { throw 'Could not determine the authenticated GitHub username.' }

$citationPath = Join-Path $root 'CITATION.cff'
$citation = Get-Content -LiteralPath $citationPath -Raw
$citation = $citation.Replace('USERNAME',$owner)
[System.IO.File]::WriteAllText($citationPath,$citation,[System.Text.UTF8Encoding]::new($false))

& (Join-Path $root 'CHECK_BEFORE_UPLOAD.ps1')
if ($LASTEXITCODE -ne 0) { throw 'Pre-upload checks failed.' }

if (-not (Test-Path -LiteralPath (Join-Path $root '.git'))) {
    git init
    git branch -M main
}

if (-not (git config user.name)) { git config user.name $owner }
if (-not (git config user.email)) { git config user.email "$owner@users.noreply.github.com" }

git add --all
$pending = git status --porcelain
if ($pending) {
    git commit -m 'Initial reproducible analysis release'
}

$remote = git remote get-url origin 2>$null
if (-not $remote) {
    $visibilityFlag = "--$Visibility"
    gh repo create $RepositoryName $visibilityFlag --description $Description --source . --remote origin --push
} else {
    git push -u origin main
}

Write-Host "Repository upload completed: $RepositoryName" -ForegroundColor Green
