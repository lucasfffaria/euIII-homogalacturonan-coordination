# Parallel MOPAC campaign runner (Windows PowerShell) used for the production run and the repeated runs
#  usage: powershell -ExecutionPolicy Bypass -File run_campaign.ps1 -Root <folder> -CampaignDir <subfolder> -MopacExe <path\to\mopac.exe>
#  - runs every <CampaignDir>\<job>\<job>.mop that has not finished yet, MaxParallel at a time
#  - automatic restart from the .arc geometry when the final gradient norm is still > GnormRestart (max 2 restarts)
#  - create a file named STOP inside the campaign folder to abort (running MOPAC processes are killed)
param(
    [int]$MaxParallel = 13,
    [string]$Root = $PSScriptRoot,
    [string]$CampaignDir = "campaign6",
    [double]$GnormRestart = 2.0,
    [int]$MaxRestarts = 2,
    [string]$MopacExe = ""
)
$ErrorActionPreference = "Continue"
$mopac = if ($MopacExe) { $MopacExe } else { Join-Path $Root "_tools\mopac-23.2.5-win\bin\mopac.exe" }   # original layout
if (-not (Test-Path $mopac)) { Write-Host "MOPAC not found at $mopac"; exit 1 }
$camp = Join-Path $Root $CampaignDir
$log = Join-Path $camp "progress.log"
$stopFile = Join-Path $camp "STOP"
if (Test-Path $stopFile) { Remove-Item $stopFile -Force }
$env:OMP_NUM_THREADS = "2"
$env:MKL_NUM_THREADS = "2"

# kill leftovers from a previous campaign
# NOTE: this kills EVERY running mopac process on the computer, including those of another runner started earlier;
# never run two instances at the same time.
Get-Process -Name mopac -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue

function Log($msg) { $line = "$(Get-Date -Format s) $msg"; Write-Host $line; Add-Content -Path $log -Value $line -Encoding ascii }

function Is-Done($outFile) {
    if (-not (Test-Path $outFile)) { return $false }
    return (Select-String -Path $outFile -Pattern "ended normally|JOB ENDED NORMALLY|== MOPAC DONE ==" -Quiet)
}
function Get-GNorm($outFile) {
    $m = Select-String -Path $outFile -Pattern "GRADIENT NORM\s+=\s+([\d.]+)" | Select-Object -Last 1
    if ($m) { return [double]$m.Matches[0].Groups[1].Value } else { return 1e9 }
}
function Make-Restart($dir, $base, $stage) {
    # build <base>_r<stage>.mop from the last .arc (block after FINAL GEOMETRY OBTAINED)
    $prev = if ($stage -eq 1) { $base } else { "$($base)_r$($stage-1)" }
    $arc = Join-Path $dir "$prev.arc"
    if (-not (Test-Path $arc)) { return $null }
    $lines = Get-Content $arc
    $idx = -1
    for ($i = 0; $i -lt $lines.Count; $i++) { if ($lines[$i] -match "FINAL GEOMETRY OBTAINED") { $idx = $i } }
    if ($idx -lt 0) { return $null }
    $block = $lines[($idx + 1)..($lines.Count - 1)]
    $newName = "$($base)_r$stage"
    $block[0] = ($block[0] -replace "^\s+", "")   # keywords line
    $block[1] = "$newName restart $stage from $prev"
    Set-Content -Path (Join-Path $dir "$newName.mop") -Value $block -Encoding ascii
    return $newName
}

$jobs = @()
Get-ChildItem -Path $camp -Directory | Sort-Object Name | ForEach-Object {
    $mop = Get-ChildItem -Path $_.FullName -Filter "*.mop" | Where-Object { $_.BaseName -notmatch "_r\d+$" } | Select-Object -First 1
    if ($mop) {
        $base = $mop.BaseName
        # determine current stage: find highest finished stage
        $stage = 0; $name = $base
        for ($s = 1; $s -le $MaxRestarts; $s++) { if (Test-Path (Join-Path $_.FullName "$($base)_r$s.mop")) { $stage = $s; $name = "$($base)_r$s" } }
        $out = Join-Path $_.FullName "$name.out"
        $done = Is-Done $out
        $needs = $true
        if ($done) {
            $g = Get-GNorm $out
            if ($g -le $GnormRestart -or $stage -ge $MaxRestarts) { $needs = $false }
            else { $stage += 1; $name = Make-Restart $_.FullName $base $stage; if (-not $name) { $needs = $false } }
        }
        if ($needs) { $jobs += [pscustomobject]@{ Dir = $_.FullName; Base = $base; Name = $name; Stage = $stage } }
    }
}
$jobs = @($jobs | Where-Object { $_.Base -like "CTRL_*" }) + @($jobs | Where-Object { $_.Base -notlike "CTRL_*" })
Log "START: $($jobs.Count) pending jobs, MaxParallel=$MaxParallel"

$running = @{}
$queue = [System.Collections.Queue]::new()
$jobs | ForEach-Object { $queue.Enqueue($_) }
$t0 = Get-Date
while ($queue.Count -gt 0 -or $running.Count -gt 0) {
    if (Test-Path $stopFile) {
        Log "STOP file found - killing $($running.Count) running jobs"
        foreach ($k in @($running.Keys)) { try { $running[$k].Proc.Kill() } catch {} }
        break
    }
    foreach ($k in @($running.Keys)) {
        $r = $running[$k]; $p = $r.Proc
        if ($p.HasExited) {
            $dt = [int]((Get-Date) - $r.Start).TotalSeconds
            $out = Join-Path $r.Job.Dir "$($r.Job.Name).out"
            $g = Get-GNorm $out
            Log "DONE  $($r.Job.Name) exit=$($p.ExitCode) gnorm=$g ${dt}s"
            $running.Remove($k)
            if ($g -gt $GnormRestart -and $r.Job.Stage -lt $MaxRestarts) {
                $ns = $r.Job.Stage + 1
                $nn = Make-Restart $r.Job.Dir $r.Job.Base $ns
                if ($nn) { $queue.Enqueue([pscustomobject]@{ Dir = $r.Job.Dir; Base = $r.Job.Base; Name = $nn; Stage = $ns }); Log "QUEUE restart $nn" }
            }
        }
    }
    while ($running.Count -lt $MaxParallel -and $queue.Count -gt 0) {
        $j = $queue.Dequeue()
        $out = Join-Path $j.Dir "$($j.Name).out"
        if (Test-Path $out) { Remove-Item $out -Force }
        $p = Start-Process -FilePath $mopac -ArgumentList "$($j.Name).mop" -WorkingDirectory $j.Dir -WindowStyle Hidden -PassThru
        $running[$j.Name] = @{ Proc = $p; Start = (Get-Date); Job = $j }
        Log "START $($j.Name) (stage $($j.Stage))"
    }
    Start-Sleep -Seconds 5
}
$el = [int]((Get-Date) - $t0).TotalMinutes
Log "ALL DONE in $el min"
