#!/usr/bin/env python3
"""Bounded sequential 3D volume showcase; preserves all other experiments."""
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/volume-showcase'
EXE=ROOT/'build/release/benchmarks/volume_showcase'

def main():
    for directory in [OUT/'raw',OUT/'visualization']:directory.mkdir(parents=True,exist_ok=True)
    cases=[(n,method,problem,.003) for problem in [1,0] for n in [6,10] for method in ['dg','ofdg']]
    cases += [(6,'ofdg',2,.003),(10,'ofdg',1,.0015)]
    env=os.environ.copy();env['PATH']='/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin'
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
    files=[ROOT/'benchmarks/volume_showcase.cpp',ROOT/'scripts/run_volume_showcase.py',EXE,*sorted((ROOT/'src').glob('*.cpp')),*sorted((ROOT/'src').glob('*.hpp'))]
    fingerprints={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    try:
        processors=int(subprocess.check_output(['/usr/sbin/sysctl','-n','hw.logicalcpu'],text=True).strip())
    except Exception: processors=None
    state={'description':'3D periodic volume advection, affine Q2 hexes; production OFDG unchanged', 'logical_processors':processors,'ranks':1,'threads_per_rank':1,'planned_cases':len(cases),'source_sha256':fingerprints,'compute_seconds':0,'peak_mib':0,'checks':[]}
    if (OUT/'run.json').exists():
        state=json.loads((OUT/'run.json').read_text())
        if state['source_sha256']!=fingerprints:raise RuntimeError('Source changed: preserve old directory and use a fresh run')
    done={c['case'] for c in state['checks']};all_rows=[]
    for check in state['checks']:
        if check['returncode']==0:all_rows.extend(csv.DictReader((OUT/'raw'/(check['case']+'.csv')).open()))
    for n,method,problem,dt in cases:
        if state['compute_seconds']>=2700:break
        key=f'n{n}-{method}-p{problem}-dt{dt}'
        if key in done:continue
        prefix=str(OUT/'visualization'/key) if problem==1 and n==10 and dt==.003 else '-'
        command=[str(EXE),str(n),method,str(problem),str(dt),'.3',prefix]
        start=time.monotonic();peak=0;failure=None
        with (OUT/'raw'/(key+'.csv')).open('w') as stdout,(OUT/'raw'/(key+'.log')).open('w') as stderr:
            process=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
            while process.poll() is None:
                entries=[tuple(map(int,line.split())) for line in subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,rss='],text=True).splitlines()]
                children={process.pid}
                while True:
                    expanded=children|{pid for pid,parent,rss in entries if parent in children}
                    if expanded==children:break
                    children=expanded
                rss=sum(rss for pid,parent,rss in entries if pid in children);peak=max(peak,rss)
                if time.monotonic()-start>300 or rss>4*1024*1024:
                    failure='time or memory ceiling';os.killpg(process.pid,signal.SIGKILL);process.wait();break
                time.sleep(.05)
        elapsed=time.monotonic()-start;state['compute_seconds']+=elapsed;state['peak_mib']=max(state['peak_mib'],peak/1024)
        state['checks'].append({'case':key,'command':command,'seconds':elapsed,'peak_mib':peak/1024,'returncode':process.returncode,'failure':failure})
        if process.returncode==0:all_rows.extend(csv.DictReader((OUT/'raw'/(key+'.csv')).open()))
        (OUT/'run.json').write_text(json.dumps(state,indent=2))
        if all_rows:
            with (OUT/'results.csv').open('w') as f:
                writer=csv.DictWriter(f,fieldnames=all_rows[0]);writer.writeheader();writer.writerows(all_rows)
        print(key,process.returncode,f'{elapsed:.2f}s',flush=True)
    print('Total seconds',state['compute_seconds'],'peak MiB',state['peak_mib'])

if __name__=='__main__':main()
