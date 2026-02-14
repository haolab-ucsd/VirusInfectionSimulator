# Pretreatment example (pt_1run)

This example demonstrates **pretreatment (priming) with a diffusible factor `dX`**
followed by a standard infection simulation.

The purpose of this example is to provide a **minimal, reproducible reference**
for how pretreatment is handled in `vsim`.

---

## What this example does

1. Generate a hexagonal lattice of cells
2. Initialize all cells with a uniform `dX` concentration
3. Run a **pretreatment phase** until `t = 50`
4. Infect the center cell
5. Continue the stochastic simulation after infection
6. Save results for **one run only** (example-scale)

This is intentionally kept to a single run to make:
- debugging easier
- output inspection faster
- logic clearer

---

## How to run

From the project root (recommended):

```bash
pip install -e .
cd examples/pt_1run
python run.py
```

## 
### Key parameters used

- Pretreatment enabled
- Pretreatment duration: t_start = 50
- Initial diffusible factor: dX0 > 0
- Number of runs: 1

All parameters actually used are saved alongside the output CSV
as a JSON file for reproducibility.

## 
### Output

The script writes CSV output files to the local out/ directory.

Example filename:
```
run000_pretreat_dX0_10_tstart_50.csv
```

Each row corresponds to a single cell at a single time point and includes:
- cell coordinates
- simulation time
- cell state
- diffusible factor levels

## 
### Notes

- This example is **not meant for performance benchmarking**
- For sweeps or multiple runs, use the main simulation utilities
- For comparison against control (no pretreatment), see wt_1run

## 
### Intended use
- Sanity check for pretreatment logic
- Template for building custom pretreatment experiments
- Reference implementation for future extensions