#!/usr/bin/env python3
"""Small local-element sweep to order ten; order-five specialization control."""
import csv,json,os
from pathlib import Path
import run_fe_mapping_derivatives as bounded
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/fe-mapping-derivatives/high-order'
def main():
 OUT.mkdir(parents=True,exist_ok=True);(OUT/'raw').mkdir(exist_ok=True)
 inputs=['benchmarks/fe_mapping_high_order.cpp','benchmarks/fe_mapping_high_order.hpp','benchmarks/fe_mapping_jets.hpp','scripts/run_fe_mapping_high_order.py','scripts/run_fe_mapping_derivatives.py','src/curved_geometry.hpp','build/release/benchmarks/fe_mapping_high_order']
 fingerprints={p:bounded.digest(ROOT/p) for p in inputs}
 state=dict(source_sha256=fingerprints,checks=[],compute_seconds=0.,peak_mib=0.,logical_processors=os.cpu_count(),ranks=1,threads_per_rank=1)
 if (OUT/'run.json').exists():
  state=json.loads((OUT/'run.json').read_text())
  if state['source_sha256']!=fingerprints:raise RuntimeError('Source changed: preserve previous high-order run first')
 done={r['case'] for r in state['checks'] if r['returncode']==0}
 env=os.environ.copy();env['PATH']='/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin'
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
 cases=[(d,n,c,10) for d in [2,3] for c in [0,1] for n in [2,8,32]]+[(d,8,1,5) for d in [2,3]]
 bounded.OUT=OUT
 for args in cases:
  key='-'.join(map(str,args))
  if key in done:continue
  total=json.loads((OUT.parent/'run.json').read_text())['compute_seconds']+state['compute_seconds']
  if total>=2700:raise RuntimeError('Combined compute budget exhausted')
  if (OUT/'raw'/(key+'.csv')).exists():raise RuntimeError('Preserve interrupted output before retry')
  rec=bounded.run_case(key,[str(ROOT/'build/release/benchmarks/fe_mapping_high_order'),*map(str,args)],env,state)
  if rec['returncode']:raise RuntimeError('Failed: '+key)
 rows=[]
 for rec in state['checks']:rows.extend(csv.DictReader((OUT/'raw'/(rec['case']+'.csv')).open()))
 with (OUT/'results.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
