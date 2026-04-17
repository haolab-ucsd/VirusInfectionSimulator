from vsim.simulation import run_many_and_save
from vsim.model import ModelParams
import json

with open("params.json") as f:
    params = ModelParams(**json.load(f))

run_many_and_save(
    params=params,
    runs=40,
    cell_iter=10000,
    t_after_infection=1200,
    n_jobs = 20,
    out_dir="out",
)
