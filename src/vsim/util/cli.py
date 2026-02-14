import argparse
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pretreat", action="store_true",
                   help="Enable pretreatment (uniform initial dX and pre-run before infection).")
    p.add_argument("--t_start", type=float, default=0.0,
                   help="Pretreatment duration (time units). Ignored if --pretreat is not set.")
    p.add_argument("--dX0", type=float, default=10.0,
                   help="Uniform initial dX value for pretreatment (default 10). Ignored if --pretreat is not set.")
    p.add_argument("--runs", type=int, default=40)
    p.add_argument("--cell_iter", type=int, default=10000)
    p.add_argument("--n_jobs", type=int, default=20)
    p.add_argument("--out_dir", type=str, default="./",
               help="Output directory for CSV files.")
    p.add_argument("--param_file", type=str, default=None,
               help="JSON file with parameter overrides (e.g., rigoe.json).")
    p.add_argument("--param_scale", type=float, default=4.1,
               help="Scale factor applied to kinetic parameters in SCALED_KEYS (default 4.1).")
    return p.parse_args()


