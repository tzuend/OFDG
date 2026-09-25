#!/usr/bin/env python3
"""Bounded sequential surface-FEM study; preserves all other experiments."""
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/sphere-surface'
EXE=ROOT/'build/release/benchmarks/sphere_surface'

def main():
    for directory in [OUT/'raw',OUT/'visualization']:directory.mkdir(parents=True,exist_ok=True)
    cases=[(family,level,p,p,0,0) for family in ['tri','quad'] for p in [1,2,3] for level in [1,2,3]]
    cases += [(family,4,2,2,0,0) for family in ['tri','quad']]
    cases += [(family,3,2,g,0,0) for family in ['tri','quad'] for g in [1,3]]
    cases += [(family,3,2,2,6,0) for family in ['tri','quad']]
    cases += [(family,2,2,2,0,1) for family in ['tri','quad']]
    env=os.environ.copy();env['PATH']='/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin'
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
    files=[ROOT/'benchmarks/sphere_surface.cpp',EXE,ROOT.parent/'mfem/examples/ex7.cpp']
    fingerprints={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    state={'description':'MFEM ex7 surface screened Poisson, u=xy on the sphere; conforming H1 FEM, no OFDG', 'logical_processors':8,'ranks':1,'threads_per_rank':1,'planned_cases':len(cases),'source_sha256':fingerprints,'compute_seconds':0,'peak_mib':0,'checks':[]}
    if (OUT/'run.json').exists():
        state=json.loads((OUT/'run.json').read_text())
        if state['source_sha256']!=fingerprints:raise RuntimeError('Source changed: preserve old directory and use a fresh run')
    done={c['case'] for c in state['checks']};all_rows=[]
    for check in state['checks']:
        if check['returncode']==0:all_rows.extend(csv.DictReader((OUT/'raw'/(check['case']+'.csv')).open()))
    for family,level,p,g,extra,constant in cases:
        if state['compute_seconds']>=2700:break
        key=f'{family}-l{level}-p{p}-g{g}-q{extra}-constant{constant}'
        if key in done:continue
        export=level==3 and p==2 and g==2 and extra==0 and constant==0
        prefix=str(OUT/'visualization'/family) if export else '-'
        command=[str(EXE),family,str(level),str(p),str(g),str(extra),str(constant),prefix]
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
