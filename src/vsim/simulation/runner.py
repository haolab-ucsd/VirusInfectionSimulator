from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional

import json
import numpy as np
import pandas as pd
import tellurium as te
from joblib import Parallel, delayed

from vsim.util.formatting import fmt_float
from vsim.model.lattice import generate_cells
from vsim.model.reaction import generate_reactions
from vsim.model.params import ModelParams
from vsim.simulation.readout import getCellInfo


def run_once(
    *,
    params: ModelParams,
    cell_iter: int,
    seed: int,
    pretreat: bool = False,
    t_start: float = 0.0,
    dX0: float = 0.0,
    t_after_infection: float = 500.0,
    n_steps: int = 100,
    n_jobs: int = 1,
    param_scale: float = 4.1,
) -> pd.DataFrame:
    """
    Run a single stochastic simulation and return a long-form dataframe
    (one row per cell per time step).
    """
    # 1) lattice
    cells = generate_cells(cell_iter)

    # 2) build antimony string (reaction network + initial conditions + params)
    S = generate_reactions(cells, dX0=dX0 if pretreat else 0.0, params=params, param_scale=param_scale)

    # 3) load into tellurium
    r = te.loada(S)
    r.integrator = "gillespie"
    r.integrator.variable_step_size = False
    r.setSeed(seed)

    # 4) optional pretreatment phase
    if pretreat and t_start > 0:
        _ = r.simulate(0, t_start, 2)  # coarse is fine; it just advances the state

    # 5) infect center then simulate
    center = cells.index([0, 0, 0])
    r[f"UI_{center}"] = 0
    r[f"I_{center}"] = 1

    sim_t0 = t_start if pretreat else 0.0
    sim_t1 = sim_t0 + t_after_infection
    result = r.simulate(sim_t0, sim_t1, n_steps)

    # 6) readout -> dataframe
    names = r.getFloatingSpeciesIds()
    result_array = np.array(result)

    out_list = []
    for t_step in range(len(result)):
        sub = Parallel(n_jobs=n_jobs)(
            delayed(getCellInfo)(c, cells, result_array, t_step, names)
            for c in range(len(cells))
        )
        out_list.extend(sub)

    df = pd.DataFrame(
        {
            "IDX": [o[0] for o in out_list],
            "X": [o[1] for o in out_list],
            "Y": [o[2] for o in out_list],
            "Z": [o[3] for o in out_list],
            "H": [o[4] for o in out_list],
            "V": [o[5] for o in out_list],
            "T": [o[6] for o in out_list],
            "T_step": [o[7] for o in out_list],
            "Df": [o[8] for o in out_list],
            "S": [o[9] for o in out_list],
            "Deathf": [o[10] for o in out_list],
        }
    )

    # time zero alignment: always shift so infection time = 0
    df["T"] = df["T"] - sim_t0

    return df


def run_many_and_save(
    *,
    params: ModelParams,
    runs: int,
    out_dir: str | Path,
    cell_iter: int,
    pretreat: bool = False,
    t_start: float = 0.0,
    dX0: float = 0.0,
    t_after_infection: float = 500.0,
    n_steps: int = 100,
    n_jobs: int = 1,
    param_scale: float = 4.1,
    save_params_each_run: bool = True,
) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for run in range(runs):
        df = run_once(
            params=params,
            cell_iter=cell_iter,
            seed=run,
            pretreat=pretreat,
            t_start=t_start,
            dX0=dX0,
            t_after_infection=t_after_infection,
            n_steps=n_steps,
            n_jobs=n_jobs,
            param_scale=param_scale,
        )

        fname = (
            f"run{run}_"
            f"{'pretreat' if pretreat else 'control'}_"
            f"dX0_{fmt_float(dX0)}_"
            f"tstart_{fmt_float(t_start)}.csv"
        )

        out_path = out_dir / fname
        df.to_csv(out_path, index=False)

        if save_params_each_run:
            used_param_path = out_dir / f"run{run}_params_used.json"
            with open(used_param_path, "w") as f:
                json.dump(asdict(params), f, indent=2, sort_keys=True)
