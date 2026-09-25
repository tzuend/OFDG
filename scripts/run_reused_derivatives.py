#!/usr/bin/env python3
"""All five methods, identical fixtures, bounded sequential fingerprinted run."""
import csv,json,os,shutil
from pathlib import Path
import run_fe_mapping_derivatives as bounded
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/reused-derivatives'
def main():
 (OUT/'raw').mkdir(parents=True,exist_ok=True);(OUT/'fixtures').mkdir(exist_ok=True)
 inputs=['benchmarks/'+s for s in ['reused_derivatives.hpp','reused_mapping_2d.cpp','reused_mapping_3d.cpp','reused_high_order.cpp','fe_mapping_jets.hpp','fe_mapping_high_order.hpp','direct_derivatives.cpp']]+['scripts/run_reused_derivatives.py','scripts/run_fe_mapping_derivatives.py','src/ofdg_core.cpp','src/curved_geometry.hpp','build/release/libofdg.a']+['build/release/benchmarks/'+s for s in ['reused_mapping_2d','reused_mapping_3d','reused_high_order']]
 hashes={p:bounded.digest(ROOT/p) for p in inputs}
 state=dict(description='Reused affine reference matrices and chain rule, with and without centering',source_sha256=hashes,checks=[],compute_seconds=0.,peak_mib=0.,ranks=1,threads_per_rank=1)
 if (OUT/'run.json').exists():
  state=json.loads((OUT/'run.json').read_text())
  if state['source_sha256']!=hashes:raise RuntimeError('Inputs changed; preserve prior run')
 else:
  for p in inputs:
   if p.startswith('build/'):continue
   dest=OUT/'source'/p;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/p,dest)
 done={r['case'] for r in state['checks'] if r['returncode']==0};env=os.environ.copy();env['PATH']='/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin'
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
 cases=[(dim+'-'+'-'.join(map(str,args)),'reused_mapping_'+dim,args) for dim,_,args in bounded.cases()]
 for d,n,c,o in [(d,n,c,10) for d in [2,3] for c in [0,1] for n in [2,8,32]]+[(d,8,1,5) for d in [2,3]]:
  key=f'high-{d}-{n}-{c}-{o}';cases.append((key,'reused_high_order',(d,n,c,o,str(OUT/'fixtures'/(key+'.txt')))))
 bounded.OUT=OUT
 for key,exe,args in cases:
  if key in done:continue
  if state['compute_seconds']>=2700:raise RuntimeError('45 minute total limit')
  if (OUT/'raw'/(key+'.csv')).exists():raise RuntimeError('Preserve incomplete output: '+key)
  rec=bounded.run_case(key,[str(ROOT/'build/release/benchmarks'/exe),*map(str,args)],env,state)
  if rec['returncode']:raise RuntimeError('Failed: '+key)
 for prefix,name in [('high-','high-order.csv'),('','results.csv')]:
  rows=[]
  for rec in state['checks']:
   if rec['case'].startswith('high-')==bool(prefix):rows.extend(csv.DictReader((OUT/'raw'/(rec['case']+'.csv')).open()))
  with (OUT/name).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
