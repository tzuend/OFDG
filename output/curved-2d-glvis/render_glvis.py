#!/usr/bin/env python3
"""Render saved MFEM snapshots using the local GLVis server (port 19916)."""
import socket, time
from pathlib import Path
root = Path(__file__).resolve().parent
for family, label, order in [('quadrilateral', 'quadrilaterals', 'Q2'), ('triangular', 'triangles', 'P2')]:
    for stage, t in [('initial', '0'), ('final', '0.5')]:
        name = f'{family}-{stage}'
        stem = root / 'data' / name
        output = root / (name + '.png')
        commands = f"""window_size 1200 1000
window_title 'Curved {label} - {stage}'
shading cool
subdivisions 8 0
keys Rjmlc
zoom 1.15
palette_name viridis
valuerange -0.02 1.02
autoscale off
plot_caption 'Curved {label} | OFDG {order} | t = {t}'
colorbar_numberformat '%.2f'
"""
        started = time.time()
        with socket.create_connection(('127.0.0.1', 19916), timeout=10) as stream:
            stream.sendall(('solution\n' + stem.with_suffix('.mesh').read_text() + '\n' + stem.with_suffix('.gf').read_text() + '\n' + commands).encode())
            time.sleep(5)
            stream.sendall(f'screenshot {output}\n'.encode())
            for _ in range(150):
                if output.exists() and output.stat().st_mtime > started:
                    break
                time.sleep(0.1)
            else:
                raise RuntimeError(f'GLVis did not save {output}')
            time.sleep(0.5)
            # Leave each rendered view open for inspection.
        print(output)
