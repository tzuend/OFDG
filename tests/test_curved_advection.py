"""Small curved-mesh convergence and partition-independence check."""
import math
from pathlib import Path
import subprocess
import sys

executable = Path(sys.argv[1] if len(sys.argv) > 1 else 'build/release/examples/advection').resolve()
for method in ('ofdg', 'ofdg-kxrcf'):
    previous = None
    serial = None
    for n in (8, 16, 32):
        arguments = [str(executable), '-d', '2', '-n', str(n), '-o', '2',
                     '-curved', '-method', method, '-tf', '.05', '-dt', '.0005', '-no-vis']
        output = subprocess.check_output(arguments, text=True, timeout=180)
        row = next(line for line in output.splitlines() if line.startswith('method='))
        values = dict(item.split('=', 1) for item in row.split())
        error = float(values['l2_error'])
        assert math.isfinite(error)
        assert float(values['conservation_drift']) < 1e-10, row
        rate = math.log2(previous / error) if previous is not None else None
        if rate is not None:
            assert rate > 2.4, (method, n, rate)
        print(method, 'n=', n, 'L2=', error, 'rate=', rate, flush=True)
        previous = error
        if n == 8:
            serial = values
    arguments[arguments.index('-n') + 1] = '8'
    output = subprocess.check_output(['mpirun', '-np', '2', *arguments], text=True, timeout=180)
    row = next(line for line in output.splitlines() if line.startswith('method='))
    parallel = dict(item.split('=', 1) for item in row.split())
    assert abs(float(parallel['l2_error']) - float(serial['l2_error'])) < 1e-10
    assert float(parallel['conservation_drift']) < 1e-10
    print(method, 'two-rank evolution agrees', flush=True)
