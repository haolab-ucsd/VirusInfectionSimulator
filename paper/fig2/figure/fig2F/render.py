import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# Data Import (dir-based condition)
# =========================
def load_runs_rig_by_dir(pattern: str):
    """
    Expected paths like:
      ../../simulation/wt/out/run0_control_dX0_0_tstart_0.csv
      ../../simulation/rigoe/out/run0_control_dX0_0_tstart_0.csv
      ../../simulation/rigko/out/run0_control_dX0_0_tstart_0.csv

    Condition is inferred from directory name: wt / rigoe / rigko
    """
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {pattern}")

    # run and tstart from filename
    # e.g. run0_control_dX0_0_tstart_0.csv
    rx = re.compile(r"run(?P<run>\d+).*_tstart_(?P<tstart>\d+)\.csv$")

    # condition from directory (supports .../wt/... , .../rigoe/... , .../rigko/...)
    # also tolerates uppercase
    cond_rx = re.compile(r"/(wt|rigoe|rigko)/", re.IGNORECASE)

    dfs = []
    for p in paths:
        m = rx.search(p)
        if not m:
            continue

        mcond = cond_rx.search(p.replace("\\", "/"))  # just in case on Windows
        if not mcond:
            continue

        run = int(m.group("run"))
        tstart = int(m.group("tstart"))

        cond_map = {"wt": "WT", "rigoe": "RIGOE", "rigko": "RIGKO"}
        cond = cond_map[mcond.group(1).lower()]

        df = pd.read_csv(p, usecols=["T", "T_step", "S"])
        df["run"] = run
        df["tstart"] = tstart
        df["cond"] = cond
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError(
            "Files matched glob but none matched expected naming pattern / condition directories."
        )
    return pd.concat(dfs, ignore_index=True)


# =========================
# 2) Mean population size calculation
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
# 3) Plot (same as before)
# =========================
def plot_wt_rigoe_rigko(
    df_mean: pd.DataFrame,
    out_pdf: str = "wt_rigoe_rigko.pdf",
    tstart: int = 0,
    show_shade: bool = False,
    use_sem: bool = False,
    show_legend: bool = True,
):
    sub = df_mean[df_mean["tstart"] == tstart].copy()
    if sub.empty:
        raise ValueError(f"No data for tstart={tstart}")

    fig, ax = plt.subplots(figsize=(4.2, 3.2), dpi=300)

    style = {
        "WT":    dict(color="black",    linestyle="-", linewidth=2.2),
        "RIGOE": dict(color="#7f7f7f",  linestyle="-", linewidth=2.2),
        "RIGKO": dict(color="#ff7f0e",  linestyle="-", linewidth=2.2),
    }
    alpha = 0.22

    for cond in ["WT", "RIGOE", "RIGKO"]:
        s = sub[sub["cond"] == cond].sort_values("T")
        if s.empty:
            continue

        T = s["T"].to_numpy()
        y = s["mean"].to_numpy()

        if use_sem:
            band = s["sd"].to_numpy() / np.sqrt(s["n"].to_numpy())
        else:
            band = s["sd"].to_numpy()

        if show_shade:
            ax.fill_between(T, y - band, y + band,
                            color=style[cond]["color"], alpha=alpha, linewidth=0)

        ax.plot(T, y, label=cond, **style[cond])

    ax.set_xlabel("Time Post Infection")
    ax.set_ylabel("Mean population size")

    x_min, x_max = sub["T"].min(), sub["T"].max()
    x_ticks = np.arange(int(np.floor(x_min/100)*100), int(np.ceil(x_max/100)*100)+1, 100)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(x) if x % 200 == 0 and x != 0 else "" for x in x_ticks])

    ax.set_ylim(0, 550)
    y_ticks = np.arange(0, 551, 50)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([str(y) if y % 100 == 0 and y != 0 else "" for y in y_ticks])

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if show_legend:
        ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved: {out_pdf}")


# =========================
# 4) 実行例
# =========================
if __name__ == "__main__":
    pattern = "../../simulation/*/out/run*_control_dX*_tstart_*.csv"
    df = load_runs_rig_by_dir(pattern)
    summ = summarize_population(df)
    mean_df = mean_over_runs(summ)

    plot_wt_rigoe_rigko(
        mean_df,
        out_pdf="FigF_repro.pdf",
        tstart=0,
        show_shade=False,
        use_sem=False,
        show_legend=False
    )
