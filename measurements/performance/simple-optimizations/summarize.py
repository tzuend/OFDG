import sys,json,csv,statistics
from pathlib import Path
root=Path('/Users/timozund/Documents/master_thesis/codex_dev');sys.path.insert(0,str(root/'scripts'))
from run_filter_performance import rows
out=root/'measurements/performance/simple-optimizations';runs=json.loads((out/'results.json').read_text());assert len(runs)==12
summary=[];data={}
for r in runs:
 ss=[s for s in rows(out/(r['name']+'.log')) if s['repeat']>=0];assert len(ss)==15
 data[(r['p'],r['curved'],r['jump'],r['version'])]=ss
 for method in ['dg','ofdg','ofdg-kxrcf']:
  samples=[s for s in ss if s['method']==method];med=statistics.median
  summary.append(dict(p=r['p'],curved=r['curved'],jump=r['jump'],version=r['version'],method=method,
    total_ms=1000*med(s['total_seconds']/s['steps'] for s in samples),
    decay_ms=1000*med(s['decay_seconds']/s['steps'] for s in samples),
    indicator_ms=1000*med(s['detector_seconds']/s['steps'] for s in samples),
    decay_share=100*med(s['decay_seconds']/s['total_seconds'] for s in samples),peak_mib=r['peak_mib']))
for key,old in data.items():
 if key[-1]!='before':continue
 new=data[(*key[:-1],'after')]
 for a in old:
  b=next(s for s in new if s['method']==a['method'] and s['repeat']==a['repeat'])
  for k in ['l2_error','mass_drift']:assert abs(a[k]-b[k])<1e-12,(key,k,a[k],b[k])
  assert abs(a['active_count']/a['steps']-b['active_count']/b['steps'])<1e-12
with (out/'summary.csv').open('w') as f:
 w=csv.DictWriter(f,summary[0]);w.writeheader();w.writerows(summary)
lines=['# Further maintainable OFDG/KXRCF optimizations','',
'This comparison starts after the earlier KXRCF face-transformation cache change. Original benchmark datasets are retained separately.', '',
'## Implemented', '',
'- **OFDG: reject inactive local faces before geometry construction.** Connectivity and the existing mask already determine whether a face will be skipped. Moving that check ahead of MFEM transformation setup avoids work without new caches, storage, or MPI control flow.',
'- **OFDG: reuse adjacent tail projections in decay.** Successive shells share S(j)u. Carry that vector to the next shell instead of recomputing it. For p≥1 this reduces projection matrix-vector products from 2p−1 to p per component per active cell: 3→2 for P2 and 5→3 for P3. Shell arithmetic and damping factors are unchanged; P0 remains unchanged.',
'- **KXRCF: precompute the static radius factor.** h_K^((p+1)/2) is now evaluated during construction. This replaces repeated powers and order lookups with a cached scalar. This is a small cleanup, not a claimed major speedup.', '',
'## Measured effect', '',
'Six configurations, 64² quadrilaterals, one rank, one library thread, five repeats plus warm-up. Identical short T=0.005 trajectories use production DG/RK4 and filtering. Old/new run order alternates. Times are median milliseconds per step; reductions compare those medians, not paired DG-normalized ratios. Small ungated differences can be ordinary timing variation. All three optimizations were measured together; no isolated speedup is claimed for the radius factor or projection reuse.', '',
'| Geometry / degree / data | OFDG phase before → after | Gated OFDG phase before → after | Full gated step before → after | Full gated reduction |',
'|---|---:|---:|---:|---:|']
for p,c,j,version in data:
 if version!='before':continue
 def get(method,version):return next(s for s in summary if (s['p'],s['curved'],s['jump'],s['method'],s['version'])==(p,c,j,method,version))
 ob,oa=get('ofdg','before'),get('ofdg','after');gb,ga=get('ofdg-kxrcf','before'),get('ofdg-kxrcf','after')
 lines.append(f"| {'Curved' if c else 'Affine'} / P{p} / {'jump' if j else 'smooth'} | {ob['decay_ms']:.2f} → {oa['decay_ms']:.2f} | {gb['decay_ms']:.2f} → {ga['decay_ms']:.2f} | {gb['total_ms']:.2f} → {ga['total_ms']:.2f} | {100*(1-ga['total_ms']/gb['total_ms']):.1f}% |")
lines+=['', '## Validation and resources', '',
'Existing OFDG fingerprints/inactive-cell checks, KXRCF tests, mixed elements, OEDG regressions, and curved geometry conservation/basis/derivative tests passed, including two-rank curved checks. Across 180 measured samples, matching old/new final L2 errors and mass drifts agree within 1e-12, and mean activation counts agree within 1e-12.', '',
f"Benchmark processes took {sum(r['wall_seconds'] for r in runs)/60:.2f} minutes; maximum sampled aggregate RSS was {max(r['peak_mib'] for r in runs):.1f} MiB. Each check stayed under five minutes and 4 GB. These changes add no growing cache storage. Raw logs, commands, binaries, source snapshots, and summary.csv are in measurements/performance/simple-optimizations.", '',
'## Opportunities deliberately not implemented', '',
f"In this new sample, the entire remaining gated OFDG phase occupies {min(s['decay_share'] for s in summary if s['version']=='after' and s['method']=='ofdg-kxrcf'):.1f}–{max(s['decay_share'] for s in summary if s['version']=='after' and s['method']=='ofdg-kxrcf'):.1f}% of total gated stepping time. Removing that entire phase is an unattainable upper bound on selective-derivative savings here; actual savings would be smaller. These proportions depend on degree, geometry, and activity.", '',
'1. **Build derivatives only where a jump needs them.** A cell must be included if it is active or neighbors an active cell; shared-face accumulation and imported neighbor states must remain consistent. This replaces the simple build-all invariant with a per-call required-cell mask and validity rules. Complexity: moderate, across derivative preparation, face traversal, and MPI tests. It could remove much derivative work for sparse activity, but cannot remove global normalization or KXRCF. Current timers combine derivative preparation with other OFDG work, so a derivative-only saving is not measured.', '',
'2. **Extend OFDG geometry caching to wave-speed evaluation and shared faces.** Local trace matrices are already cached, but wave-speed evaluation still obtains quadrature rules and evaluates normals; shared faces also rebuild trace data and neighbor derivative states. A unified cache could avoid these operations. Complexity: moderate: reconcile the existing trace cache with speed quadrature (including the OEDG rule distinction), own stable transformations, handle imported elements and face orientation, and account for additional memory. It could save a material part of face-processing time, but that part is not separately timed. The earlier KXRCF cache speedup is evidence that geometry setup can matter, not a prediction of the OFDG saving.', '',
'3. **Tensor-product evaluation for quadrilaterals/hexahedra.** Replace dense high-order evaluation with successive one-dimensional contractions. Potential is greatest at higher order; P1–P3 timings do not establish a worthwhile gain. Complexity: high: specialized tensor kernels, basis/order handling, and a separate dense path for simplices and mixed meshes. This would be a larger numerical-library change, not a small readability-neutral optimization.', '',
'4. **Reduce temporary allocations in shared face physics.** Prescribed-advection evaluation constructs two velocity vectors at each quadrature point; Euler constructs primitive-state vectors. Caller-owned reusable scratch could avoid these allocations. Complexity: moderate and cross-cutting: extend the FacePhysics interface and its implementations/call sites, preserving ownership and concurrent-use semantics. A separate transport-only interface could also avoid unused spectral-speed work for KXRCF but adds another physics contract to maintain. Allocation and spectral-speed costs are not separately profiled, so no percentage saving is justified. This is a candidate for targeted profiling if KXRCF becomes a bottleneck.', '',
'No reduction of curved quadrature, freezing of state-dependent inflow, or changed filtering cadence was made. Those would alter numerical behaviour rather than merely eliminate redundant work. With the simple fixes in place, further optimization should start with targeted profiling and a demonstrated application need.']
(root/'docs/further-optimizations.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
