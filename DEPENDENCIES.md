# Reproducible dependency baseline

The validated dependency baseline for development, CI, and future numerical
measurements is:

- MFEM commit `a1ce49fb5742ad3c46e3e3732e5f9b04e4ca6e18`;
- MFEM configured with `MFEM_USE_MPI=YES MFEM_USE_METIS_5=YES`;
- C++17;
- Open MPI, Hypre, and METIS;
- Python 3 with NumPy and Matplotlib;
- `latexmk` for report builds;
- Ubuntu 24.04 for the continuous-integration reference environment.

The exact MFEM commit is deliberate. CI must fetch that object directly rather
than building the moving `master` branch. The study runner additionally records
the project revision, MFEM revision, compiler, host, and command in every raw
row. Package-manager patch versions may differ between local machines, so final
measurements must retain the recorded metadata and must not mix hosts within a
timing comparison.

The existing Make interface remains authoritative. A compatible local checkout
is expected at `../mfem`, or may be selected with:

```sh
make MFEM_DIR=/absolute/path/to/mfem test-release
```

Before producing an approved dataset, verify both revisions:

```sh
git rev-parse HEAD
git -C ../mfem rev-parse HEAD
```
