"""Focused preview checks: physical export, exact data, limits, and MPI export."""
import csv,json,math,sys,tempfile,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_preview import ROOT,bounded,expand
manifest=json.loads((ROOT/'experiments/preview_manifest.json').read_text())
assert len(expand(manifest,8))==100
assert max(r['ranks'] for r in expand(manifest,8))==2
assert max(r['ranks'] for r in expand(manifest,1))==1
with tempfile.TemporaryDirectory() as directory:
    d=Path(directory)
    r=bounded([sys.executable,'-c','import time; time.sleep(2)'],d/'timeout.log',.25,4096)
    assert r['status']=='time_limit',r
    r=bounded([sys.executable,'-c','import time; x=bytearray(30_000_000); time.sleep(2)'],d/'memory.log',5,10)
    assert r['status']=='memory_limit',r
    subprocess.run([str(ROOT/'build/release/tests/test_preview')],cwd=d,check=True)
    rows=list(csv.DictReader((d/'preview_export_check.rank000000.csv').open()))
    assert len([r for r in rows if r['sample']=='polynomial'])==21
    assert abs(float(next(r for r in rows if r['sample']=='cell_average')['value'])-1/3)<1e-13
    # T=0 tests the new piecewise initial data away from its two jump cells.
    cmd=[str(ROOT/'build/release/examples/advection'),'-piecewise','-n','20','-o','3','-tf','0','-method','dg','-no-vis']
    subprocess.run(cmd+['-profile',str(d/'initial')],check=True,stdout=subprocess.DEVNULL)
    rows=list(csv.DictReader((d/'initial.rank000000.csv').open()))
    def exact(x):
        x=x-math.floor(x)
        return math.sin(2*math.pi*x) if .3<=x<=.8 else math.cos(2*math.pi*x)-.5
    for row in rows:
        x=float(row['x'])
        if row['sample']=='polynomial' and min(abs(x-.3),abs(x-.8),abs(x),abs(x-1))>.055:
            assert abs(float(row['value'])-exact(x))<1e-5,(x,row['value'],exact(x))
    # Actual exported means must agree when the same curved test is partitioned.
    cmd=[str(ROOT/'build/release/examples/advection'),'-d','2','-curved','-n','8','-o','2','-tf','.02','-method','ofdg-kxrcf','-no-vis']
    for ranks in [1,2]:
        subprocess.run((['mpirun','-np','2'] if ranks==2 else [])+cmd+['-profile',str(d/f'np{ranks}')],check=True,stdout=subprocess.DEVNULL)
    def means(prefix):
        return sorted((float(r['x']),float(r['y']),float(r['value'])) for p in d.glob(prefix+'*.csv') for r in csv.DictReader(p.open()) if r['sample']=='cell_average')
    a,b=means('np1'),means('np2')
    assert len(a)==len(b)==64
    assert max(abs(x-y) for left,right in zip(a,b) for x,y in zip(left,right))<1e-11
print('Preview runner, physical exports, piecewise data, and serial/MPI checks passed')
