"""
Minimal example:
- WT condition
- no pretreatment
- single run
- write one CSV
"""

from pathlib import Path

from vsim.model import ModelParams
from vsim.simulation import run_many_and_save


def main():
    out_dir = Path("./output")
    out_dir.mkdir(exist_ok=True)

    # WT default parameters
    params = ModelParams()

    # run exactly once with minimal settings
    run_many_and_save(
        params=params,
        runs=1,
        out_dir=out_dir,
        cell_iter=1000,      # small for quick test
        pretreat=False,
        t_start=0.0,
        dX0=0.0,
        n_jobs=1,
        param_scale=4.1,
    )


if __name__ == "__main__":
    main()
