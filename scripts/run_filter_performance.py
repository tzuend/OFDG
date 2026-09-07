#!/usr/bin/env python3
"""Bounded single-rank performance experiment, separate from numerical studies."""
import argparse,csv,hashlib,json,subprocess,time
from pathlib import Path
from run_preview import bounded,available_cpus,ENV
ROOT=Path(__file__).resolve().parents[1]
BIN=ROOT/'build/release/benchmarks/filter_performance'

def rows(log):
    data=[]
    for line in log.read_text().splitlines():
        if not line.startswith('BENCH '):continue
        values=dict(s.split('=',1) for s in line.split()[1:])
        data.append({k:(v if k=='method' else float(v)) for k,v in values.items()})
    return data

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--trajectory-time',type=float,default=.005)
    ap.add_argument('--prior-budget-seconds',type=float,default=300.)
    ap.add_argument('--output',type=Path,default=ROOT/'measurements/performance/2026-09-07')
    a=ap.parse_args(); out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    (out/'raw').mkdir(exist_ok=True)
    paths=[BIN,ROOT/'build/release/libofdg.a',Path(__file__),ROOT/'benchmarks/filter_performance.cpp',ROOT.parent/'mfem/libmfem.a']
    fingerprint={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    configs=[dict(n=n,p=2,curved=g,jump=False) for n in [16,32,64,128] for g in [False,True]]
    configs += [dict(n=n,p=p,curved=g,jump=False) for p in [1,3] for n in [16,32,64] for g in [False,True]]
    configs += [dict(n=n,p=2,curved=g,jump=True) for n in [16,64] for g in [False,True]]
    f=out/'results.json'
    state=json.loads(f.read_text()) if f.exists() else dict(fingerprint=fingerprint,cpus=available_cpus(),ranks=1,
        thread_limits={k:ENV[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']},
        elapsed_seconds=a.prior_budget_seconds,prior_budget_seconds=a.prior_budget_seconds,trajectory_time=a.trajectory_time,runs=[],validation=[],configs=configs,limits=dict(memory_mib=4096,check_seconds=300,launch_seconds=2700,total_seconds=3600),
        project_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        mfem_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT.parent/'mfem',text=True).strip())
    assert state['fingerprint']==fingerprint,'Use a new directory after changing benchmark sources or binaries.'
    if state.pop('pending',None):state['elapsed_seconds']+=300
    def save():
        f.write_text(json.dumps(state,indent=2)+'\n')
        samples=[{**{k:run[k] for k in ['id','status','peak_mib']},**row} for run in state['runs'] for row in run.get('samples',[])]
        if samples:
            with (out/'samples.csv').open('w') as s:
                writer=csv.DictWriter(s,list(samples[0]));writer.writeheader();writer.writerows(samples)
    def execute(command,identifier,group):
        state['pending']=identifier;save()
        log=out/'raw'/f'{identifier}.log'
        result=bounded(command,log,min(300,3600-state['elapsed_seconds']),4096)
        state.pop('pending');state['elapsed_seconds']+=result['wall_seconds']
        result.update(id=identifier,command=command);state[group].append(result);save()
        return result,log
    # Verify that the instrumented path agrees with the ordinary advection driver.
    if not state.get('validated'):
        for curved in [False,True]:
            flag='-curved' if curved else '-affine'
            result,log=execute([str(BIN),'-n','8','-o','2',flag,'-r','1','-ms','0'],'verify_bench_'+str(curved),'validation')
            assert result['status']=='ok',result
            measured={r['method']:r for r in rows(log) if r['repeat']==0}
            for method in ['dg','ofdg','ofdg-kxrcf']:
                cmd=[str(ROOT/'build/release/examples/advection'),'-d','2','-n','8','-o','2','-tf','.05','-c','.15','-s','4','-method',method,'-no-vis']
                if curved:cmd+=['-curved']
                ordinary,_=execute(cmd,f'verify_driver_{curved}_{method}','validation')
                assert ordinary['status']=='ok',ordinary
                assert abs(float(ordinary['metrics']['l2_error'])-measured[method]['l2_error'])<1e-10,(ordinary,measured)
        state['validated']=True;save()
        print('Driver equivalence checks passed for both geometries and all three methods.',flush=True)
    for cfg in configs:
        identifier=f"{'curved' if cfg['curved'] else 'affine'}_p{cfg['p']}_n{cfg['n']}_{'jump' if cfg['jump'] else 'smooth'}"
        if any(r['id']==identifier for r in state['runs']):continue
        if state['elapsed_seconds']>=2700:break
        command=[str(BIN),'-n',str(cfg['n']),'-o',str(cfg['p']),'-curved' if cfg['curved'] else '-affine','-jump' if cfg['jump'] else '-smooth','-tf',str(a.trajectory_time)]
        print('Running '+identifier,flush=True)
        result,log=execute(command,identifier,'runs');result['samples']=rows(log);save()
        print(f"  {result['status']}: {result['wall_seconds']:.1f}s, sampled peak {result['peak_mib']:.0f} MiB; total {state['elapsed_seconds']/60:.1f} min",flush=True)
    save()
if __name__=='__main__':main()
