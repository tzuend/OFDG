#!/usr/bin/env python3
"""Sequential, checkpointed laptop study. No changes to full-study assets.

Every solver, reference and adaptive rerun consumes the same persisted budget.
Resume requires identical sources/configuration; use a new output directory after edits.
"""
import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from run_study import parse_metrics

ROOT = Path(__file__).resolve().parents[1]
ENV = {**os.environ, 'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1',
       'VECLIB_MAXIMUM_THREADS':'1', 'MKL_NUM_THREADS':'1'}


def available_cpus():
    try:
        return int(subprocess.check_output(['/usr/sbin/sysctl','-n','hw.logicalcpu'], text=True))
    except (OSError, ValueError, subprocess.CalledProcessError):
        return 1


def tree_rss(pid):
    text = subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,rss='], text=True)
    entries = [tuple(map(int, line.split())) for line in text.splitlines()]
    descendants = {pid}
    while True:
        expanded = descendants | {p for p, parent, _ in entries if parent in descendants}
        if expanded == descendants:
            return sum(rss for p, _, rss in entries if p in descendants) / 1024
        descendants = expanded


def bounded(command, log, seconds, memory_mib):
    start = time.monotonic()
    peak = 0.
    reason = ''
    with log.open('w') as output:
        process = subprocess.Popen(command, cwd=log.parent, env=ENV,
            stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while process.poll() is None:
                peak = max(peak, tree_rss(process.pid))
                if peak > memory_mib:
                    reason = 'memory_limit'
                elif time.monotonic() - start >= seconds:
                    reason = 'time_limit'
                if reason:
                    os.killpg(process.pid, signal.SIGKILL)
                    break
                time.sleep(.2)
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
            process.wait()
    elapsed = time.monotonic() - start
    return dict(status=reason or ('ok' if process.returncode == 0 else 'failed'),
                exit_code=process.returncode, wall_seconds=elapsed, peak_mib=peak,
                metrics=parse_metrics(log.read_text()))


def expand(manifest, cpus):
    runs=[]
    for case in manifest['cases']:
        pairs = zip(case['orders'],case['resolutions']) if case.get('zip_order_resolution') else itertools.product(case['orders'],case['resolutions'])
        for order,n in pairs:
            methods=case['methods'] + (['oedg'] if n==case.get('oedg_resolution') else [])
            for method in methods:
                ranks = min(2, max(1,cpus//2)) if case.get('dimension',1)==2 and n>=16 else 1
                runs.append(dict(case=case['name'], program=case['program'],order=order,n=n,
                                 method=method,ranks=ranks,cfl=case['cfl'],tf=case['final_time'],args=case['args']))
    return runs


def run_id(run):
    return f"{run['case']}_{run['method']}_p{run['order']}_n{run['n']}_c{run['cfl']:g}"


def source_fingerprint():
    paths=[]
    for directory in ['src','examples','benchmarks']:
        paths += sorted((ROOT/directory).rglob('*.cpp')) + sorted((ROOT/directory).rglob('*.hpp'))
    paths += [ROOT/'Makefile',Path(__file__),ROOT/'experiments/preview_manifest.json']
    digest=hashlib.sha256()
    for p in paths:
        digest.update(str(p.relative_to(ROOT)).encode()); digest.update(p.read_bytes())
    for program in ['examples/advection','examples/burgers','examples/euler','benchmarks/euler_reference_weno']:
        digest.update((ROOT/'build/release'/program).read_bytes())
    return digest.hexdigest()


def save(path, state):
    temp=path.with_suffix('.tmp'); temp.write_text(json.dumps(state,indent=2)+'\n'); temp.replace(path)
    rows=[]
    for row in state['runs']:
        rows.append({**{k:v for k,v in row.items() if k not in ['metrics','args','command']},**row.get('metrics',{})})
    if rows:
        with path.with_name('results.csv').open('w') as stream:
            writer=csv.DictWriter(stream,sorted(set().union(*(r.keys() for r in rows))))
            writer.writeheader(); writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'measurements/study/preview')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'experiments/preview_manifest.json').read_text())
    cpus=available_cpus(); runs=expand(manifest,cpus)
    if args.dry_run:
        print(json.dumps(dict(cpus=cpus,runs=runs),indent=2)); return 0
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    for name in ['raw','profiles','references']:(output/name).mkdir(exist_ok=True)
    statefile=output/'results.json'
    fingerprint=source_fingerprint()
    state=json.loads(statefile.read_text()) if statefile.exists() else dict(fingerprint=fingerprint,cpus=cpus,elapsed_seconds=0.,runs=[],references=[])
    if state['fingerprint']!=fingerprint:
        raise RuntimeError('Sources or binaries changed: choose a new output directory.')
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (output/'source.patch').write_text(subprocess.check_output(['git','diff','HEAD'],cwd=ROOT,text=True))
    state['provenance']={name:subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip() for name,path in [('project',ROOT),('mfem',ROOT.parent/'mfem')]}
    state['thread_limits']={k:ENV[k] for k in ENV if k.endswith('NUM_THREADS') or k=='VECLIB_MAXIMUM_THREADS'}
    limits=manifest['limits']
    # Charge an interrupted check conservatively before resuming its unfinished run.
    if state.pop('pending', None) is not None:
        state['elapsed_seconds'] += limits['check_seconds']
        save(statefile,state)
    def execute(command, stem):
        remaining=limits['total_seconds']-state['elapsed_seconds']
        if state['elapsed_seconds']>=limits['launch_seconds']:
            return dict(status='skipped_budget',wall_seconds=0.,peak_mib=0.,metrics={})
        state['pending']=stem
        save(statefile,state)
        result=bounded(command,output/'raw'/f'{stem}.log',min(limits['check_seconds'],remaining),limits['memory_mib'])
        state.pop('pending',None)
        state['elapsed_seconds']+=result['wall_seconds']
        return result
    done={r['id'] for r in state['runs']}
    def solve(run):
        stem=run_id(run)
        if stem in done:return
        prefix=output/'profiles'/stem
        command=[str(ROOT/'build/release/examples'/run['program']),*map(str,run['args']),'-n',str(run['n']),'-o',str(run['order']),'-tf',str(run['tf']),'-c',str(run['cfl']),'-s','4','-method',run['method'],'-cadence','auto','-no-vis','-profile',str(prefix)]
        if run['program']=='euler':command+=['-r','0','-positivity','-no-save']
        if run['ranks']>1:command=['mpirun','-np',str(run['ranks']),*command]
        print(stem,flush=True)
        if any(r['case']==run['case'] and r['method']==run['method'] and
               r['n']<run['n'] and r['status'] in ('time_limit','memory_limit') for r in state['runs']):
            result=dict(status='skipped_larger_after_limit',wall_seconds=0.,peak_mib=0.,metrics={})
        else:
            result=execute(command,stem)
        state['runs'].append({**run,'id':stem,'command':command,'profile_prefix':str(prefix),**result})
        done.add(stem);save(statefile,state)
        print(f"  {result['status']} {result['wall_seconds']:.2f}s; total {state['elapsed_seconds']/60:.1f}min",flush=True)
    # Scalar cases precede references; references precede Euler comparisons.
    for run in runs:
        if run['program']=='euler':break
        solve(run)
    for reference in manifest['references']:
        if any(r['name']==reference['name'] for r in state['references']):continue
        dest=output/'references'/f"{reference['name']}.csv"
        command=[str(ROOT/'build/release/benchmarks/euler_reference_weno'),'--problem',str(reference['problem']),'--cells',str(reference['cells']),'--final-time',str(reference['final_time']),'--output',str(dest)]
        print('reference '+reference['name'],flush=True)
        state['references'].append({**reference,'command':command,**execute(command,reference['name'])})
        save(statefile,state)
    for run in runs:solve(run)
    # Temporal sensitivity uses the same final time and exact-error quadrature.
    for method in manifest['cases'][0]['methods']:
        base=next(r for r in runs if r['case']=='smooth_advection' and r['order']==3 and r['n']==256 and r['method']==method)
        half={**base,'cfl':base['cfl']/2}
        solve(half)
        original=next(r for r in state['runs'] if r['id']==run_id(base))
        refined=next(r for r in state['runs'] if r['id']==run_id(half))
        if original['status']==refined['status']=='ok':
            a=float(original['metrics']['l2_error']);b=float(refined['metrics']['l2_error'])
            change=abs(a-b)/max(abs(b),1e-300)
            state.setdefault('temporal_sensitivity',{})[method]=change
            if change>.05:
                for r in runs:
                    if r['case']=='smooth_advection' and r['order']==3 and r['method']==method:
                        solve({**r,'cfl':r['cfl']/2})
    save(statefile,state)
    print(f"Completed preview: {state['elapsed_seconds']/60:.2f} minutes; {output}")
    return 0

if __name__=='__main__':sys.exit(main())
