# Perturbation example (pert_1run)

This example demonstrates a **single-run perturbation simulation**
in which model parameters are modified via an external JSON file.

It is intended as a **template for parameter perturbation experiments**
(e.g. knockouts, pathway inhibition, sensitivity analysis).

---

## What this example does

1. Generate a hexagonal lattice of cells
2. Load parameter overrides from a JSON file
3. Apply those overrides to the default WT parameter set
4. Infect the center cell at `t = 0`
5. Run a stochastic infection simulation
6. Save results for **one run only** (example-scale)

In this example, we reproduce the **IFNR knockout–like condition**
reported in the paper by disabling X-mediated effects.

---

## How to run

From the project root (recommended):

```bash
pip install -e .
cd examples/pert_1run
python run.py
```
## 
### Parameter specification
Perturbations are specified via a JSON file using vsim parameter names.

Example: ifnrko.json
```
{
  "k_dXI": 0.0,
  "k_Xmax": 0.0
}
```

These correspond to disabling:
- X-mediated infection attenuation
- X-dependent death signaling

All unspecified parameters remain at their WT default values.
##
### Key parameters used
- No pretreatment
- Infection starts at t = 0
- Initial diffusible factor: dX0 = 0
- Number of runs: 1
- Parameters overridden via JSON

The exact parameters used are saved as ``params_used.json`` in the output directory for reproducibility.
## 
### Output

The script writes CSV output files to the local out/ directory.

Example filename:
```
ifnrko_run000_control_dX0_0_tstart_0.csv
```

Each row corresponds to a single cell at a single time point and includes:
- cell coordinates
- simulation time
- cell state
- diffusible factor levels
## 
### Notes
- This example shows how to decouple simulation logic from parameter choice
- JSON-based perturbations make simulations reproducible and scriptable
- The same structure can be used for knockouts, partial inhibition, or scans
## 
### Intended use
- Parameter perturbation experiments
- Reproducing knockout-like conditions (e.g. IFNRKO)
- Template for more systematic perturbation studies