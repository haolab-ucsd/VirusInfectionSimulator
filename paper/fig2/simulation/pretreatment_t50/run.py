from vsim.simulation import run_many_and_save
from vsim.model import ModelParams
import json

params = ModelParams()
run_many_and_save(
    params=params,
    runs=40,
    cell_iter=10000,
    pretreat = True,
    t_start = 50.0,
    t_after_infection=1200,
    dX0 = 1.0,
    n_jobs = 20,
    out_dir="out",
)
