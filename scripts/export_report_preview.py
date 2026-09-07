#!/usr/bin/env python3
"""Export thesis tables from retained preview measurements; never launch solvers."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'measurements/study/preview'
DEST = ROOT / 'report/generated/preview'
LABEL = {'dg': 'DG', 'ofdg': 'OFDG', 'ofdg-kxrcf': 'OFDG--KXRCF', 'oedg': 'OEDG'}


def table(name, spec, header, rows, caption, label):
    text = ('% Generated from retained preview data. Do not edit by hand.\n'
            '\\begin{table}[htbp]\n\\centering\\small\n'
            '\\begin{tabular}{@{}' + spec + '@{}}\n\\toprule\n' +
            ' & '.join(header) + '\\\\\\midrule\n')
    text += ''.join(' & '.join(row) + '\\\\\n' for row in rows)
    text += ('\\bottomrule\n\\end{tabular}\n\\caption{' + caption + '}\n'
             '\\label{' + label + '}\\end{table}\n')
    (DEST / (name + '.tex')).write_text(text)


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    a = json.loads((SOURCE / 'analysis.json').read_text())
    state = json.loads((SOURCE / 'results.json').read_text())
    report_analysis = json.loads((DEST / 'analysis.json').read_text())
    if report_analysis != a:
        raise ValueError('Report analysis differs from retained preview analysis; reconcile the dataset before exporting claims')
    rows = []
    for case, title, degrees in [('smooth_advection', 'Advection', [1,2,3]),
                                 ('smooth_burgers', 'Burgers', [2]),
                                 ('straight_advection_2d', 'Straight 2D', [2]),
                                 ('curved_advection_2d', 'Curved 2D', [2])]:
        for degree in degrees:
            for method in list(LABEL)[:3]:
                series = [v for v in a['convergence'] if v['case']==case and v['p']==degree and v['method']==method]
                v = max(series, key=lambda v:v['n'])
                rows.append([title + f' ({degree})', str(v['n']), LABEL[method],
                             f"{v['l2']:.3e}", f"{v['observed_order']:.2f}"])
    table('smooth_table', 'll lrr', ['Problem (degree)', '$N$', 'Method', '$L^2$ error', 'Rate'], rows,
          'Finest-grid smooth errors and rates on the last refinement. For 2D, $N$ is cells per direction and degree 2 denotes $Q^2$. Rates above the nominal order are finite-range observations.', 'tab:smooth-errors')
    rows = [[('Transport' if v['case']=='piecewise_advection' else 'Burgers'), LABEL[v['method']],
             f"{v['minimum']:.5f}", f"{v['maximum']:.5f}"] for v in a['extrema']]
    table('extrema_table', 'llrr', ['Problem', 'Method', 'Minimum', 'Maximum'], rows,
          'Sampled polynomial extrema on the 256-cell scalar cases. Initial ranges are $[-1,\\sin(0.6\\pi)]$ for transport and $[-0.5,1.5]$ for Burgers. Values rounded to five decimals can hide small excursions.', 'tab:scalar-extrema')
    rows = [[{'lax':'Lax','shu_osher':'Shu--Osher'}[case], f"{v['mean_absolute_density_difference']:.6f}", f"{v['max_density_difference']:.6f}"] for case,v in a['reference_sensitivity'].items()]
    table('reference_table', 'lrr', ['Problem', 'Mean absolute difference', 'Maximum difference'], rows,
          'Density sensitivity of 2,048- versus 4,096-cell WENO references after conservative restriction to the coarse grid. These are differences between references, not DG errors.', 'tab:reference-sensitivity')
    runs = [r for r in state['runs'] if r['case']=='riemann_2d']
    rows = []
    for n in [32,64]:
        for method in list(LABEL)[:3]:
            v = next(r for r in runs if r['n']==n and r['method']==method)['metrics']
            rows.append([str(n), LABEL[method], f"{float(v['min_density']):.3f}", f"{float(v['min_pressure']):.3f}", v['limited_elements']])
    table('positivity_table', 'rlrrr', ['$N$', 'Method', 'Min. density', 'Min. pressure', 'Limited events'], rows,
          'Reported final sampled minima and accumulated cell-limiting events for 2D Riemann. All rows have zero rejected steps. Positive sampled values do not certify positivity at every point.', 'tab:riemann-positivity')
    inputs = [SOURCE/'results.json', SOURCE/'analysis.json', ROOT/'experiments/preview_manifest.json', ROOT/'scripts/analyze_preview.py', Path(__file__)]
    inputs += sorted((SOURCE/'profiles').glob('*.csv')) + sorted((SOURCE/'references').glob('*.csv'))
    provenance = {'purpose':'Preliminary thesis figures and tables; no new simulations',
                  'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (DEST/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')

if __name__ == '__main__':
    main()
