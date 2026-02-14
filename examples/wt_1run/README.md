# Control example (wt_1run)

This example demonstrates a **baseline infection simulation without pretreatment**.

It serves as the **control (WT) reference** against which pretreatment
or parameter-modified simulations can be compared.

---

## What this example does

1. Generate a hexagonal lattice of cells
2. Initialize all cells with **no pretreatment** (`dX = 0`)
3. Infect the center cell at `t = 0`
4. Run a stochastic infection simulation
5. Save results for **one run only** (example-scale)

This example represents the **default WT condition** used throughout the project.

---

## How to run

From the project root (recommended):

```bash
pip install -e .
cd examples/wt_1run
python run.py
```
## 
### Key parameters used
- No pretreatment
- Infection starts at t = 0
- Initial diffusible factor: dX0 = 0
- Number of runs: 1

All parameters actually used are saved alongside the output CSV
as a JSON file for reproducibility.

## 
### Output
The script writes CSV output files to the local out/ directory.

Example filename:
```
run000_control_dX0_0_tstart_0.csv
```

Each row corresponds to a single cell at a single time point and includes:
- cell coordinates
- simulation time
- cell state
- diffusible factor levels

## 
### Notes
- This example is intended as a baseline reference
- Output from this example is directly comparable to pt_1run
- For parameter sweeps or multiple stochastic runs, use the main simulation APIs
## 
### Intended use
- Baseline (WT) reference simulation
- Comparison against pretreatment or knockout conditions
- Minimal working example for infection dynamics