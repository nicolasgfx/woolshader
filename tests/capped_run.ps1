# capped_run.ps1 -- run a command under a hard wall-clock cap (seamless plan M0).
#
#   powershell -NoProfile -ExecutionPolicy Bypass -File tests\capped_run.ps1 -CapSeconds 3600 -Log run.log `
#       -Exe build\knit_check.exe -ArgString "--scene wool-swatch"
#
# The command's stdout goes to -Log, stderr to <Log>.err, the cap record to <Log>.cap. If it is
# still running at the cap, its whole process tree is killed (taskkill /T /F), "FAIL(time cap)
# after N s" is written to the record and the exit code is 124. Otherwise the command's exit code
# is returned. The cooperative caps inside knit_check (knit.time_cap_s, --part-cap) end a run
# cleanly at their checkpoints; this script is the hard limit for whatever does not reach a
# checkpoint in time (a long device call, a library operation, a lock wait).
param(
    [Parameter(Mandatory = $true)][double]$CapSeconds,
    [Parameter(Mandatory = $true)][string]$Log,
    [Parameter(Mandatory = $true)][string]$Exe,
    [string]$ArgString = ""
)
$start = Get-Date
"capped_run: cap $CapSeconds s; started $($start.ToString('s')); command: $Exe $ArgString" | Out-File -FilePath "$Log.cap" -Encoding utf8
$sp = @{ FilePath = $Exe; NoNewWindow = $true; PassThru = $true; RedirectStandardOutput = $Log; RedirectStandardError = "$Log.err" }
if ($ArgString -ne "") { $sp.ArgumentList = $ArgString }
$p = Start-Process @sp
$null = $p.Handle  # keep the handle so ExitCode stays readable
$done = $p.WaitForExit([int]([Math]::Min($CapSeconds * 1000, [int]::MaxValue)))
$elapsed = ((Get-Date) - $start).TotalSeconds
if (-not $done) {
    & taskkill.exe /T /F /PID $p.Id 2>&1 | Out-Null
    $p.WaitForExit(10000) | Out-Null
    $line = "FAIL(time cap) after {0:N1} s (cap {1} s): process tree of PID {2} killed" -f $elapsed, $CapSeconds, $p.Id
    $line | Out-File -FilePath "$Log.cap" -Append -Encoding utf8
    Write-Output $line
    exit 124
}
$code = $p.ExitCode
("finished after {0:N1} s, exit code {1}" -f $elapsed, $code) | Out-File -FilePath "$Log.cap" -Append -Encoding utf8
exit $code
