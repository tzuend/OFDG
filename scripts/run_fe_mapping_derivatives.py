#!/usr/bin/env python3
"""Bounded, fingerprinted FE mapping experiment; no solver changes."""
import argparse,csv,hashlib,json,os,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/fe-mapping-derivatives'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cases():
 c=[(g,4,p,0.,'gll','coupled',0) for g in ['quad','tri'] for p in [2,3,4]]
 c += [(g,n,p,a,'gll',m,0) for m,a in [('separable',.7),('coupled',1.2)] for g in ['quad','tri'] for p in [2,3,4] for n in [4,8]]
 c += [(g,4,3,a,b,m,0) for g in ['quad','tri'] for m,a in [('coupled',1.2),('separable',.7)] for b in ['gl','bernstein']]
 c += [(g,4,3,a,'gll',m,8) for g in ['quad','tri'] for m,a in [('coupled',1.2),('separable',.7)]]
 c += [(g,8,3,a,'gll',m,0) for g in ['quad','tri'] for m,a in [('coupled',.4),('separable',.3)]]
 c += [(g,16,3,a,'gll',m,0) for g in ['quad','tri'] for m,a in [('coupled',1.2),('separable',.7)]]
 for args in c:yield '2d','fe_mapping_derivatives',args
 for curved in [0,1]:
  for n in [2,4]:
   for p in [2,3]:yield '3d','fe_mapping_derivatives_3d',(n,p,curved,0)
 yield '3d','fe_mapping_derivatives_3d',(4,3,1,8)
def run_case(key,command,env,state):
 raw=OUT/'raw';start=time.monotonic();peak=0;reason=None
 with (raw/(key+'.csv')).open('w') as stdout,(raw/(key+'.log')).open('w') as stderr:
  proc=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
  while proc.poll() is None:
   entries=[tuple(map(int,s.split())) for s in subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,rss='],text=True).splitlines()];children={proc.pid}
   while True:
    found=children|{pid for pid,parent,rss in entries if parent in children}
    if found==children:break
    children=found
   rss=sum(rss for pid,parent,rss in entries if pid in children);peak=max(peak,rss)
   if time.monotonic()-start>300 or rss>4*1024*1024:
    reason='time or memory ceiling';os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
   time.sleep(.1)
 elapsed=time.monotonic()-start
 state['compute_seconds']+=elapsed;state['peak_mib']=max(state['peak_mib'],peak/1024)
 record=dict(case=key,command=command,seconds=elapsed,peak_mib=peak/1024,returncode=proc.returncode,failure=reason)
 state['checks'].append(record);(OUT/'run.json').write_text(json.dumps(state,indent=2))
 print(f'{key}: status={proc.returncode}, {elapsed:.2f}s, {peak/1024:.1f} MiB',flush=True)
 return record
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--dimension',choices=['2d','3d','all'],default='all');args=parser.parse_args()
 (OUT/'raw').mkdir(parents=True,exist_ok=True)
 inputs=['benchmarks/fe_mapping_jets.hpp','benchmarks/fe_mapping_derivatives.cpp','benchmarks/fe_mapping_derivatives_3d.cpp','benchmarks/direct_derivatives.cpp','src/curved_geometry.hpp','scripts/run_fe_mapping_derivatives.py','build/release/benchmarks/fe_mapping_derivatives','build/release/benchmarks/fe_mapping_derivatives_3d']
 fingerprints={p:digest(ROOT/p) for p in inputs}
 state=dict(description='Three-method FE geometry derivative comparison',source_sha256=fingerprints,checks=[],compute_seconds=0.,peak_mib=0.,logical_processors=os.cpu_count(),ranks=1,threads_per_rank=1,mfem_config=(ROOT.parent/'mfem/config/config.mk').read_text())
 if (OUT/'run.json').exists():
  state=json.loads((OUT/'run.json').read_text())
  if state['source_sha256']!=fingerprints:raise RuntimeError('Inputs changed: preserve previous results before rerunning')
 done={r['case'] for r in state['checks'] if r['returncode']==0}
 env=os.environ.copy();env['PATH']='/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin'
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
 for dim,exe,values in cases():
  if args.dimension!='all' and args.dimension!=dim:continue
  key=dim+'-'+'-'.join(map(str,values))
  if key in done:continue
  if state['compute_seconds']>=2700:raise RuntimeError('Total compute ceiling reached')
  if (OUT/'raw'/(key+'.csv')).exists():raise RuntimeError('Preserve failed/interrupted output before retry: '+key)
  record=run_case(key,[str(ROOT/'build/release/benchmarks'/exe),*map(str,values)],env,state)
  if record['returncode']:raise RuntimeError('Case failed: '+key)
 rows=[]
 for c in state['checks']:
  if c['returncode']==0:rows.extend(csv.DictReader((OUT/'raw'/(c['case']+'.csv')).open()))
 if rows:
  with (OUT/'results.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 print('Successful configurations:',len(done)+sum(1 for r in state['checks'] if r['case'] not in done),flush=True)
if __name__=='__main__':main()
