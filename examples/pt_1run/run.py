from pathlib import Path

from vsim.model import ModelParams
from vsim.simulation.runner import run_many_and_save

def main() -> None:
    out_dir = Path("./out")
    out_dir.mkdir(parents=True, exist_ok=True)

    params = ModelParams()

    run_many_and_save(
        params=params,
        out_dir=out_dir,
        runs=1,
        cell_iter=10000,
        pretreat=True,           # Pretreatment ON
        t_start=50.0,
        dX0=10.0,                # Initial dX in pretreatment
        t_after_infection=500.0,
        n_steps=100,
        n_jobs=20,
        param_scale=4.1,
    )

if __name__ == "__main__":
    main()
