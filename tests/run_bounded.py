"""Run a check with a five-minute/2-GiB process-tree ceiling (POSIX)."""
import os
import signal
import subprocess
import sys
import time

start = time.monotonic()
process = subprocess.Popen(sys.argv[1:], start_new_session=True)
peak = 0
try:
    while process.poll() is None:
        rows = subprocess.check_output(['ps', '-axo', 'pid=,ppid=,rss='], text=True)
        entries = [tuple(map(int, line.split())) for line in rows.splitlines()]
        descendants = {process.pid}
        while True:
            expanded = descendants | {pid for pid, parent, _ in entries if parent in descendants}
            if expanded == descendants:
                break
            descendants = expanded
        rss = sum(memory for pid, _, memory in entries if pid in descendants)
        peak = max(peak, rss)
        if time.monotonic() - start > 300 or rss > 2 * 1024 * 1024:
            raise RuntimeError('Check exceeded five minutes or 2 GiB resident memory')
        time.sleep(0.2)
except BaseException:
    os.killpg(process.pid, signal.SIGKILL)
    process.wait()
    raise
print(f'Bounded check: {time.monotonic()-start:.2f}s, peak tree RSS {peak/1024:.1f} MiB', flush=True)
sys.exit(process.returncode)
