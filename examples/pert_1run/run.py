from __future__ import annotations

import json
from pathlib import Path

from vsim.model import ModelParams
from vsim.util.params_io import load_params_json, override_params
from vsim.simulation.runner import run_many_and_save


def main() -> None:
    here = Path(__file__).resolve().parent
    out_dir = here / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    param_path = here / "ifnrko.json"
    updates = load_params_json(str(param_path))

    base = ModelParams()
    params = override_params(base, updates)

    with open(out_dir / "params_used.json", "w") as f:
        json.dump(params.__dict__, f, indent=2, sort_keys=True)

    run_many_and_save(
        params=params,
        out_dir=out_dir,
        runs=1,
        cell_iter=10000,
        pretreat=False,
        t_start=0.0,
        dX0=0.0,
        t_after_infection=500.0,
        n_steps=100,
        n_jobs=20,
        param_scale=4.1
    )

    print(f"[pert_1run] Done. Outputs are in: {out_dir}")


if __name__ == "__main__":
    main()
