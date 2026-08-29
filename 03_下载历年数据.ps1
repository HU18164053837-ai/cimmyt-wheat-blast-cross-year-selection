$ErrorActionPreference = 'Stop'
$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$historyRoot = Join-Path $scriptRoot 'history'
$rawRoot = Join-Path $historyRoot 'raw_downloads'
$metaRoot = Join-Path $historyRoot 'metadata'
New-Item -ItemType Directory -Force -Path $rawRoot,$metaRoot | Out-Null

$datasets = @(
    [pscustomobject]@{Key='cycle_2023'; Handle='11529/10549099'; Scope='2023 cycle: 14HLBSN, 14HZAN, 40SAWSN, 55IBWSN plus companion nurseries'},
    [pscustomobject]@{Key='cycle_2022'; Handle='11529/10548927'; Scope='2022 cycle: 13HLBSN, 13HZAN, 38SAWSN, 53IBWSN, 54IBWSN plus companion nurseries'},
    [pscustomobject]@{Key='HLBSN_2018_2021'; Handle='11529/10548688'; Scope='8th to 12th HLBSN'},
    [pscustomobject]@{Key='IBWSN_2018_2021'; Handle='11529/10548690'; Scope='50th to 52nd IBWSN'},
    [pscustomobject]@{Key='SAWSN_2018_2021'; Handle='11529/10548691'; Scope='35th to 37th SAWSN'}
)

$manifest = @()
foreach ($ds in $datasets) {
    $datasetDir = Join-Path $rawRoot $ds.Key
    New-Item -ItemType Directory -Force -Path $datasetDir | Out-Null
    $landing = 'https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:' + $ds.Handle
    $page = Invoke-WebRequest -Uri $landing -TimeoutSec 60 -UseBasicParsing
    $ids = [regex]::Matches($page.Content, 'datafile/\d+') | ForEach-Object { $_.Value.Split('/')[-1] } | Sort-Object -Unique
    if (-not $ids) { throw "No file IDs discovered for $($ds.Handle)" }
    $existingFiles = @(Get-ChildItem -LiteralPath $datasetDir -File | Where-Object { $_.Extension -ne '.tmp' })
    if ($existingFiles.Count -ge $ids.Count) {
        foreach ($f in $existingFiles) {
            $manifest += [pscustomobject]@{
                dataset_key=$ds.Key; handle=$ds.Handle; scope=$ds.Scope; datafile_id='resolved_previous_run'; file_name=$f.Name;
                bytes=$f.Length; md5=(Get-FileHash -LiteralPath $f.FullName -Algorithm MD5).Hash.ToLowerInvariant();
                download_url=$landing; downloaded_at=(Get-Date).ToString('s')
            }
        }
        continue
    }
    foreach ($id in $ids) {
        $tmp = Join-Path $datasetDir ("download_" + $id + '.tmp')
        $downloadUrl = 'https://data.cimmyt.org/api/access/datafile/' + $id
        $response = $null
        for ($attempt=1; $attempt -le 4; $attempt++) {
            try {
                $response = Invoke-WebRequest -Uri $downloadUrl -TimeoutSec 120 -OutFile $tmp -PassThru -UseBasicParsing
                break
            } catch {
                if ($attempt -eq 4) { throw }
                Start-Sleep -Seconds ([Math]::Pow(2,$attempt))
            }
        }
        $cd = [string]$response.Headers['Content-Disposition']
        if ($cd -match "filename\*=UTF-8''([^;]+)") {
            $fileName = [uri]::UnescapeDataString($Matches[1].Trim('"'))
        } elseif ($cd -match 'filename="?([^";]+)"?') {
            $fileName = $Matches[1]
        } else {
            $fileName = 'datafile_' + $id + '.bin'
        }
        $dest = Join-Path $datasetDir $fileName
        if (Test-Path -LiteralPath $dest) { Remove-Item -LiteralPath $tmp -Force } else { Move-Item -LiteralPath $tmp -Destination $dest }
        $f = Get-Item -LiteralPath $dest
        $manifest += [pscustomobject]@{
            dataset_key=$ds.Key; handle=$ds.Handle; scope=$ds.Scope; datafile_id=$id; file_name=$f.Name;
            bytes=$f.Length; md5=(Get-FileHash -LiteralPath $f.FullName -Algorithm MD5).Hash.ToLowerInvariant();
            download_url=$downloadUrl; downloaded_at=(Get-Date).ToString('s')
        }
    }
}
$manifest | Sort-Object dataset_key,file_name | Export-Csv -LiteralPath (Join-Path $metaRoot 'download_manifest.csv') -NoTypeInformation -Encoding UTF8
$datasets | Export-Csv -LiteralPath (Join-Path $metaRoot 'target_datasets.csv') -NoTypeInformation -Encoding UTF8
Write-Output ("Downloaded {0} files from {1} datasets." -f $manifest.Count,$datasets.Count)
