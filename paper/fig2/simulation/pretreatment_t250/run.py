from vsim.simulation import run_many_and_save
from vsim.model import ModelParams
import json

params = ModelParams()
run_many_and_save(
    params=params,
    runs=40,
    cell_iter=10000,
    pretreat = True,
    t_start = 250.0,
    dX0 = 10.0,
    n_jobs = 20,
    out_dir="out",
)
