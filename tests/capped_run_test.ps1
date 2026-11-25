# CTest capped_run_kills_and_reports: a process (and its child) that outlives a 2 s cap is killed
# within seconds, the run exits 124 and the record says FAIL(time cap); a quick command keeps
# its own exit code.
param([Parameter(Mandatory = $true)][string]$Script)
$tmp = Join-Path ([IO.Path]::GetTempPath()) ("capped_run_test_" + [Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $tmp | Out-Null
$log = Join-Path $tmp "run.log"
$t0 = Get-Date
# A parent that starts a long-sleeping child: the whole tree must go.
& powershell -NoProfile -ExecutionPolicy Bypass -File $Script -CapSeconds 2 -Log $log -Exe "cmd.exe" -ArgString "/c start /b ping -n 60 127.0.0.1 >nul & ping -n 60 127.0.0.1 >nul"
$code = $LASTEXITCODE
$elapsed = ((Get-Date) - $t0).TotalSeconds
$cap = Get-Content "$log.cap" -Raw
if ($code -ne 124) { Write-Output "FAIL: exit code $code, expected 124"; exit 1 }
if ($elapsed -gt 30) { Write-Output "FAIL: took $elapsed s"; exit 1 }
if ($cap -notmatch 'FAIL\(time cap\)') { Write-Output "FAIL: no time-cap record: $cap"; exit 1 }
& powershell -NoProfile -ExecutionPolicy Bypass -File $Script -CapSeconds 30 -Log $log -Exe "cmd.exe" -ArgString "/c exit 7"
if ($LASTEXITCODE -ne 7) { Write-Output "FAIL: quick command exit $LASTEXITCODE, expected 7"; exit 1 }
Remove-Item -Recurse -Force $tmp -ErrorAction SilentlyContinue
Write-Output ("[ OK ] capped_run_kills_and_reports (killed after {0:N1} s)" -f $elapsed)
exit 0
