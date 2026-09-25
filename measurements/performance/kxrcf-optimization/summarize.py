import sys,json,statistics,csv
from pathlib import Path
root=Path('/Users/timozund/Documents/master_thesis/codex_dev')
sys.path.insert(0,str(root/'scripts'))
from run_filter_performance import rows
out=root/'measurements/performance/kxrcf-optimization'
runs=json.loads((out/'results.json').read_text()); summary=[]; indexed={}
assert len(runs)==12
for r in runs:
    samples=[s for s in rows(out/(r['name']+'.log')) if s['repeat']>=0]
    indexed[(r['n'],r['p'],r['curved'],r['jump'],r['version'])]=samples
    base={s['repeat']:s for s in samples if s['method']=='dg'}
    gated=[s for s in samples if s['method']=='ofdg-kxrcf']
    med=statistics.median
    summary.append(dict(n=r['n'],p=r['p'],curved=r['curved'],jump=r['jump'],version=r['version'],
        detector_ms=1000*med(s['detector_seconds']/s['steps'] for s in gated),
        total_ms=1000*med(s['total_seconds']/s['steps'] for s in gated),
        slowdown=med((s['total_seconds']/s['steps'])/(base[s['repeat']]['total_seconds']/base[s['repeat']]['steps']) for s in gated),
        detector_share=100*med(s['detector_seconds']/s['total_seconds'] for s in gated),peak_mib=r['peak_mib']))
for key,old in indexed.items():
    if key[-1]!='before':continue
    new=indexed[(*key[:-1],'after')]
    for a in old:
        b=next(s for s in new if s['method']==a['method'] and s['repeat']==a['repeat'])
        assert abs(a['l2_error']-b['l2_error'])<1e-12
        assert abs(a['mass_drift']-b['mass_drift'])<1e-12
        assert abs(a['active_count']/a['steps']-b['active_count']/b['steps'])<1e-12
with (out/'summary.csv').open('w') as f:
    w=csv.DictWriter(f,summary[0]);w.writeheader();w.writerows(summary)
lines=['# KXRCF optimization check','',
'KXRCF computes no solution derivatives. The earlier combined filtering share included both the indicator and the remaining OFDG work. The indicator itself still evaluates cell norms and face traces at quadrature points and determines inflow. Curved overintegration increases those costs.', '',
'The implementation previously requested MFEM face transformations on every indicator call. That reconstructs element and face geometry from the mesh. The change retains owned transformations alongside the existing static face cache, using MFEM’s supported overload. Integration points, solution traces, and velocity/physics are still evaluated on each call. Public interfaces, quadrature, thresholds, and arithmetic of the indicator are unchanged. Reconstruct the indicator after mesh changes, as before. Shared MPI faces retain their existing path.', '',
'## Before/after measurements','',
'One rank, one library thread, five measured repetitions plus warm-up; unchanged RK4 benchmark at T=0.005. The order of old/new binaries alternates between configurations. Reported times are medians per step. Total speed changes include run-to-run processor/cache variation; direct indicator timing is the more targeted measurement. Original results remain separate.', '',
'| Geometry / degree / grid / data | KX before → after (ms) | KX reduction | Total gated before → after (ms) | Gated / DG before → after |',
'|---|---:|---:|---:|---:|']
for a in summary:
    if a['version']!='before':continue
    b=next(v for v in summary if v['version']=='after' and all(v[k]==a[k] for k in ['n','p','curved','jump']))
    label=f"{'Curved' if a['curved'] else 'Affine'} / P{a['p']} / {a['n']}² / {'jump' if a['jump'] else 'smooth'}"
    lines.append(f"| {label} | {a['detector_ms']:.2f} → {b['detector_ms']:.2f} | {100*(1-b['detector_ms']/a['detector_ms']):.1f}% | {a['total_ms']:.2f} → {b['total_ms']:.2f} | {a['slowdown']:.3f}× → {b['slowdown']:.3f}× |")
lines+=['', '## Checks and limits','',
'All 180 measured samples agree with their matching old/new counterpart in final L2 error and mass drift within 1e-12, and average activation counts within 1e-12. Existing KXRCF fingerprints, mixed-element tests, and curved geometry tests in serial and two ranks pass. A dedicated regression checks spatially varying, changing velocity and reuse after unrelated mesh transformation lookups.', '',
f"Numerical benchmark elapsed: {sum(r['wall_seconds'] for r in runs)/60:.2f} minutes. Maximum sampled aggregate RSS: {max(r['peak_mib'] for r in runs):.1f} MiB. Every process completed within 300 seconds and 4 GB. The CSV includes peak RSS for each old/new run: these are whole-process peaks, not an isolated cache allocation measurement.", '',
'The tradeoff is additional per-face geometry storage. No change to OFDG derivative preparation was made. Avoiding derivatives on inactive cells and their unneeded neighbors is a possible next optimization, but requires careful local/MPI dependency handling. Reducing curved quadrature or caching state-dependent inflow would affect numerical behaviour and is not part of this change.', '',
'Raw logs, exact commands, source snapshots, executable hashes, and timing CSV are in measurements/performance/kxrcf-optimization. The original benchmark brief describes the pre-optimization implementation; this document records the subsequent change.']
(root/'docs/kxrcf-optimization.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
