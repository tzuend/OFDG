# Curved 2D advection — GLVis images

Four native GLVis PNGs (2400 × 1690), showing initial t=0 and final t=0.5 states on curved quadrilateral and triangular meshes. No report files were edited.

## Problem and numerical settings

Periodic unit square; scalar advection velocity (0.7, 0.3); initial data sin²(π(x+y)). Exact transported solution sin²(π(x+y−t)). We chose t=0.5 to show a visible shift: at t=1 this wave returns to its initial state.

Start with 12 × 12 squares: 144 quadrilaterals or 288 triangles. Geometry degree 3, formed by nodally interpolating the map (x,y) ↦ (x+d,y+d), d=0.04 sin(2πx) sin(2πy). This fixes the outer square boundary. Solution degree 2: Q2 on quadrilaterals (1296 DOFs), P2 on triangles (1728 DOFs). Physical L² initialization, RK4, CFL 0.15, OFDG after each full timestep, 153 steps. These are ungated OFDG runs; the method selection overrides the legacy use-kxrcf option printed in the option listing.

Both runs were sequential, one MPI rank, one numerical-library thread, monitored with the existing five-minute/2-GiB process-tree guard. Measured total bounded runtime: about 3 seconds for both simulations. Solver times: 0.437 s (quads), 1.198 s (triangles). L² errors: 0.00929237 and 0.00704832, respectively. See raw logs under data/.

## How they were run

An isolated copy of examples/advection/example_advection.cpp is saved at tmp/curved-visualizations/advection_frames.cpp. It only adds full-precision initial/final mesh and solution exports (and adjusts the include path). Production solver sources were not changed.

From the codex_dev directory:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make -f Makefile -f output/curved-2d-glvis/frames.mk build/release/advection_frames
python3 tests/run_bounded.py build/release/advection_frames -d 2 -n 12 -o 2 -s 4 -c 0.15 -tf 0.5 -method ofdg -curved -quadrilateral -no-vis -profile output/curved-2d-glvis/data/quadrilateral
python3 tests/run_bounded.py build/release/advection_frames -d 2 -n 12 -o 2 -s 4 -c 0.15 -tf 0.5 -method ofdg -curved -triangular -no-vis -profile output/curved-2d-glvis/data/triangular
```

No mpirun is necessary for these single-process runs.

## How they were rendered

Sent each saved .mesh/.gf pair as a solution stream to the existing GLVis 4.5 server at localhost:19916. GLVis rendered its own PNG screenshots. Settings: top-down orthographic view, lighting off, element edges on, curved/DG shading (cool), subdivision level 8, viridis palette, shared colour range [−0.02,1.02]. The colour bar represents the scalar solution u; the horizontal/vertical image directions are x/y on [0,1]. No smoothing across DG element boundaries is imposed by the export.

To render again while GLVis is listening on port 19916:

```sh
python3 output/curved-2d-glvis/render_glvis.py
```

The script leaves the generated GLVis windows open. Raw snapshots allow re-rendering without rerunning the simulations. This is a visual example, not a convergence study or a comparison between stabilization methods.
