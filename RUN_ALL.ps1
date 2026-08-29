param(
  [string]$Python = "python",
  [string]$Rscript = "Rscript"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = $PSScriptRoot
$env:CIMMYT_ANALYSIS_ROOT = $ProjectRoot
$LogDir = Join-Path $ProjectRoot "运行记录"
New-Item -ItemType Directory -Force $LogDir | Out-Null
$LogFile = Join-Path $LogDir "end_to_end_rerun.log"
$Manifest = Join-Path $LogDir "run_manifest.csv"
$Checksums = Join-Path $LogDir "key_file_checksums.csv"

"Start: $(Get-Date -Format o)" | Set-Content $LogFile -Encoding UTF8
"Project root: $ProjectRoot" | Add-Content $LogFile -Encoding UTF8

$steps = @()
$steps += [pscustomobject]@{Order=4; Type="Python"; File=(Get-ChildItem $ProjectRoot -Filter "04_*.py" | Select-Object -First 1).FullName}
foreach($n in 5..14){
  $steps += [pscustomobject]@{Order=$n; Type="R"; File=(Get-ChildItem $ProjectRoot -Filter (("{0:D2}_*.R" -f $n)) | Select-Object -First 1).FullName}
}

$records = @()
foreach($step in $steps){
  if(-not $step.File){ throw "Missing script for step $($step.Order)" }
  $start = Get-Date
  "RUN $($step.Order): $([IO.Path]::GetFileName($step.File))" | Tee-Object -FilePath $LogFile -Append
  if($step.Type -eq "Python"){
    & $Python $step.File *>> $LogFile
  } else {
    & $Rscript $step.File *>> $LogFile
  }
  $exit = $LASTEXITCODE
  $records += [pscustomobject]@{
    order=$step.Order
    script=[IO.Path]::GetFileName($step.File)
    sha256=(Get-FileHash $step.File -Algorithm SHA256).Hash
    exit_status=$exit
    elapsed_seconds=[math]::Round(((Get-Date)-$start).TotalSeconds,3)
  }
  if($exit -ne 0){
    $records | Export-Csv $Manifest -NoTypeInformation -Encoding UTF8
    throw "Step $($step.Order) failed with exit status $exit. See $LogFile"
  }
}

$summary = Import-Csv (Join-Path $ProjectRoot "history\cross_year_analysis\tables\00_analysis_summary.csv")
$portfolio = Import-Csv (Join-Path $ProjectRoot "history\validation_portfolio_analysis\tables\00_analysis_summary.csv")
"KEY COUNTS" | Add-Content $LogFile -Encoding UTF8
$summary | ForEach-Object { "$($_.metric)=$($_.value)" | Add-Content $LogFile -Encoding UTF8 }
$portfolio | ForEach-Object { "$($_.metric)=$($_.value)" | Add-Content $LogFile -Encoding UTF8 }
"End: $(Get-Date -Format o)" | Add-Content $LogFile -Encoding UTF8
$records | Export-Csv $Manifest -NoTypeInformation -Encoding UTF8
$keyFiles = @(
  "history\standardized\historical_long_primary_observed.csv",
  "history\cross_year_analysis\tables\00_analysis_summary.csv",
  "history\validation_portfolio_analysis\tables\01_nonredundant_validation_portfolio.csv",
  "history\candidate_sensitivity_baseline_analysis\tables\07_full_funnel_multiverse.csv",
  "history\candidate_sensitivity_baseline_analysis\tables\08_full_funnel_candidate_stability.csv",
  "history\temporal_uncertainty_negative_control\tables\00_analysis_summary.csv",
  "history\submission_consistency_audit\tables\submission_consistency_audit.csv"
)
$keyFiles | ForEach-Object {
  $full = Join-Path $ProjectRoot $_
  if(Test-Path -LiteralPath $full){
    [pscustomobject]@{relative_path=$_; bytes=(Get-Item -LiteralPath $full).Length; sha256=(Get-FileHash -LiteralPath $full -Algorithm SHA256).Hash}
  }
} | Export-Csv $Checksums -NoTypeInformation -Encoding UTF8
Write-Output "End-to-end run completed. Log: $LogFile"
