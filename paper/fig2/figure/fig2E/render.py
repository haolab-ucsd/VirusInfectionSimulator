import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =========================
# 1) Import (cond from directory)
# =========================
def load_runs_from_dirs(control_pattern: str, ifnrko_pattern: str):
    # Ex:
    # ../../simulation/wt/out/run0_control_dX0_0_tstart_0.csv
    # ../../simulation/ifnrko/out/run0_control_dX0_0_tstart_0.csv
    paths_control = sorted(glob.glob(control_pattern))
    paths_ifnrko  = sorted(glob.glob(ifnrko_pattern))

    if not paths_control and not paths_ifnrko:
        raise FileNotFoundError(
            "No files matched either pattern:\n"
            f"  control: {control_pattern}\n"
            f"  ifnrko : {ifnrko_pattern}"
        )

    # Extract run and tstart from filename
    # run0_control_dX0_0_tstart_0.csv
    rx = re.compile(r"^run(?P<run>\d+)_.*?_tstart_(?P<tstart>\d+)\.csv$")

    dfs = []

    def _load_one(p: str, cond: str):
        fname = Path(p).name
        m = rx.search(fname)
        if not m:
            return  # skip silently (or raise if you prefer)

        run = int(m.group("run"))
        tstart = int(m.group("tstart"))

        df = pd.read_csv(p, usecols=["T", "T_step", "S"])
        df["run"] = run
        df["tstart"] = tstart
        df["cond"] = cond
        dfs.append(df)

    for p in paths_control:
        _load_one(p, "Control")
    for p in paths_ifnrko:
        _load_one(p, "IFNRKO")

    if not dfs:
        raise FileNotFoundError(
            "Files matched glob but none matched filename regex:\n"
            f"  control: {control_pattern}\n"
            f"  ifnrko : {ifnrko_pattern}"
        )

    return pd.concat(dfs, ignore_index=True)

# =========================
# 2) Mean population size calculation (unchanged)
# =========================
def summarize_population(df: pd.DataFrame):
    g = df.groupby(["cond", "tstart", "run", "T_step", "T"])["S"]
    out = g.agg(population=lambda s: np.sum(s.isin([1, 2, 4]))).reset_index()
    return out

def mean_over_runs(summ: pd.DataFrame):
    grp = summ.groupby(["cond", "tstart", "T_step", "T"])
    mean = grp["population"].mean()
    sd   = grp["population"].std(ddof=1)
    n    = grp.size()

    out = pd.concat([mean, sd, n], axis=1).reset_index()
    out.columns = ["cond", "tstart", "T_step", "T", "mean", "sd", "n"]
    return out

# =========================
# 3) Plot (unchanged)
# =========================
def plot_control_vs_ifnrko(
    df_mean: pd.DataFrame,
    out_pdf: str = "control_vs_ifnrko.pdf",
    tstart: int = 0,
    use_sem: bool = False,
    show_shade: bool = True,
    show_legend: bool = True,
):
    sub = df_mean[df_mean["tstart"] == tstart].copy()
    if sub.empty:
        raise ValueError(f"No data for tstart={tstart}")

    fig, ax = plt.subplots(figsize=(4.2, 3.2), dpi=300)

    style = {
        "Control": dict(color="black", linestyle="-",  linewidth=2),
        "IFNRKO":  dict(color="black", linestyle="--", linewidth=2),
    }

    alpha = 0.22

    for cond, st in style.items():
        s = sub[sub["cond"] == cond].sort_values("T")
        if s.empty:
            continue

        T = s["T"].to_numpy()
        y = s["mean"].to_numpy()

        band = (s["sd"].to_numpy() / np.sqrt(s["n"].to_numpy())) if use_sem else s["sd"].to_numpy()

        if show_shade:
            ax.fill_between(T, y - band, y + band, color=st["color"], alpha=alpha, linewidth=0)

        ax.plot(T, y, label=cond, **st)

    ax.set_xlabel("Time")
    ax.set_ylabel("Mean population size")
    ax.set_xlim(sub["T"].min(), sub["T"].max())
    ax.set_ylim(0, 550)

    y_ticks = np.arange(0, 551, 50)
    ax.set_yticks(y_ticks)
    y_labels = [str(y) if y % 100 == 0 and y != 0 else "" for y in y_ticks]
    ax.set_yticklabels(y_labels)

    x_min, x_max = sub["T"].min(), sub["T"].max()
    x_ticks = np.arange(int(np.floor(x_min / 100) * 100), int(np.ceil(x_max / 100) * 100) + 1, 100)
    ax.set_xticks(x_ticks)
    x_labels = [str(x) if x % 200 == 0 and x != 0 else "" for x in x_ticks]
    ax.set_xticklabels(x_labels)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if show_legend:
        ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved: {out_pdf}")

# =========================
# 4) Execution
# =========================
if __name__ == "__main__":
    control_pattern = "../../simulation/wt/out/run*_control_*_tstart_*.csv"
    ifnrko_pattern  = "../../simulation/ifnrko/out/run*_control_*_tstart_*.csv"

    df = load_runs_from_dirs(control_pattern, ifnrko_pattern)
    summ = summarize_population(df)
    mean_df = mean_over_runs(summ)

    plot_control_vs_ifnrko(
        mean_df,
        out_pdf="fig2E_repro.pdf",
        tstart=0,
        use_sem=False,
        show_shade=False,
        show_legend=False
    )
