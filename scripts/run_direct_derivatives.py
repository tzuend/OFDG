#!/usr/bin/env python3
"""Sequential, bounded operator-only study. No PDE solves or production edits."""
import csv, hashlib, io, json, os, signal, subprocess, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/direct-derivatives'
EXE=ROOT/'build/release/benchmarks/direct_derivatives'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'raw').mkdir(exist_ok=True)
    previous=json.loads((OUT/'run.json').read_text()) if (OUT/'run.json').exists() else None
    env=os.environ.copy();env['PATH']='/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin'
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[name]='1'
    cases=[(g,4,p,0.,'gll','coupled',0) for g in ['quad','tri'] for p in [2,3,4]]
    cases += [(g,n,p,a,'gll',m,0) for m,a in [('separable',.7),('coupled',1.2)] for g in ['quad','tri'] for p in [2,3,4] for n in [4,8]]
    cases += [(g,4,3,a,b,m,0) for g in ['quad','tri'] for m,a in [('coupled',1.2),('separable',.7)] for b in ['gl','bernstein']]
    cases += [(g,4,3,a,'gll',m,8) for g in ['quad','tri'] for m,a in [('coupled',1.2),('separable',.7)]]
    cases += [(g,8,3,a,'gll',m,0) for g in ['quad','tri'] for m,a in [('coupled',.4),('separable',.3)]]
    cases += [(g,16,3,a,'gll',m,0) for g in ['quad','tri'] for m,a in [('coupled',1.2),('separable',.7)]]
    state={'description':'Matched DG input: unchanged production recursion versus direct physical differentiation with inverse Taylor jets',
           'ranks':1,'threads_per_rank':1,'cpu_detection':'8 logical processors confirmed; one rank/thread chosen for small operator checks',
           'checks':[],'source_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [ROOT/'benchmarks/direct_derivatives.cpp',ROOT/'src/curved_geometry.hpp',EXE]},
           'mfem_config':(ROOT.parent/'mfem/config/config.mk').read_text(),'compute_seconds':0,'derivative_evaluation':'degree-three inverse Taylor jets of original mapped DG function','planned_cases':len(cases)}
    if previous: state=previous
    all_rows=[];peak_all=state.get('peak_mib',0)*1024
    done={c['case'] for c in state['checks']}
    for c in state['checks']:
        if c['returncode']==0:all_rows+=list(csv.DictReader((OUT/'raw'/f"{c['case']}.csv").open()))
    for g,n,p,a,b,m,extra in cases:
        if state['compute_seconds']>=2700:break
        key=f'{m}-{g}-n{n}-p{p}-a{a}-{b}-extra{extra}';command=[str(EXE),g,str(n),str(p),str(a),b,m,str(extra)]
        if key in done:continue
        if (OUT/'raw'/f'{key}.csv').exists():
            (OUT/'raw'/f'{key}.csv').rename(OUT/'raw'/f'{key}.interrupted.csv')
            (OUT/'raw'/f'{key}.log').rename(OUT/'raw'/f'{key}.interrupted.log')
        start=time.monotonic();peak=0;reason=None
        with (OUT/'raw'/f'{key}.csv').open('w') as stdout,(OUT/'raw'/f'{key}.log').open('w') as stderr:
            proc=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
            while proc.poll() is None:
                data=subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,rss='],text=True)
                entries=[tuple(map(int,s.split())) for s in data.splitlines()];children={proc.pid}
                while True:
                    found=children|{pid for pid,parent,rss in entries if parent in children}
                    if found==children:break
                    children=found
                rss=sum(rss for pid,parent,rss in entries if pid in children);peak=max(peak,rss)
                if time.monotonic()-start>300 or rss>4*1024*1024:
                    reason='time or memory limit';os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
                time.sleep(.1)
        elapsed=time.monotonic()-start;state['compute_seconds']+=elapsed;peak_all=max(peak_all,peak)
        state['checks'].append({'case':key,'command':command,'seconds':elapsed,'peak_mib':peak/1024,'returncode':proc.returncode,'failure':reason})
        state['peak_mib']=peak_all/1024
        if proc.returncode==0:all_rows+=list(csv.DictReader((OUT/'raw'/f'{key}.csv').open()))
        (OUT/'run.json').write_text(json.dumps(state,indent=2))
        with (OUT/'results.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=all_rows[0].keys());writer.writeheader();writer.writerows(all_rows)
        print(f'{key}: status={proc.returncode}, {elapsed:.2f}s, {peak/1024:.1f} MiB',flush=True)
    if all_rows:
        with (OUT/'results.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=all_rows[0].keys());writer.writeheader();writer.writerows(all_rows)
    print(f'Total {state["compute_seconds"]:.2f}s; {len(all_rows)} rows; peak {peak_all/1024:.1f} MiB')
if __name__=='__main__':main()
