import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_runs(pattern):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {pattern}")

    dfs = []
    for p in paths:
        m = re.search(r"run(\d+)_tstart(\d+)\.csv", p)
        run = int(m.group(1)) if m else None

        df = pd.read_csv(p, usecols=["T", "T_step", "S"])
        df["run"] = run if run is not None else p
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)

def summarize_counts(df):
    # cumulative infection = I + DA + DE (S in {1,2,4})
    g = df.groupby(["run", "T_step", "T"])["S"]

    out = g.agg(
        infected_death=lambda s: np.sum(s == 2),
        altruistic_death=lambda s: np.sum(s == 3),
        cumulative_infection=lambda s: np.sum(s.isin([1, 2, 4])),
    ).reset_index()

    return out

def mean_and_band(summ):
    grp = summ.groupby(["T_step", "T"])

    mean = grp[["infected_death", "altruistic_death", "cumulative_infection"]].mean()
    sd   = grp[["infected_death", "altruistic_death", "cumulative_infection"]].std(ddof=1)
    n    = grp.size()

    mean = mean.reset_index()
    sd   = sd.reset_index()
    n    = n.reset_index(name="n")

    merged = mean.merge(sd, on=["T_step", "T"], suffixes=("_mean", "_sd")).merge(n, on=["T_step","T"])
    return merged

def plot_fig2c(merged, out_pdf="fig2c_repro.pdf", use_sem=False):
    T = merged["T"].to_numpy()
    n = merged["n"].to_numpy()

    def band(y_mean, y_sd, n):
        if use_sem:
            y_band = y_sd / np.sqrt(n)
        else:
            y_band = y_sd
        return y_mean - y_band, y_mean + y_band

    # ===== Colour Definition =====
    COL_INFECTED = "#d62728"   # 赤
    COL_ALTRUIST = "#2ca02c"   # 緑
    COL_CUMUL    = "#000000"   # 黒

    SHADE_ALPHA = 0.25

    # ===== Data  =====
    red_m = merged["infected_death_mean"].to_numpy()
    red_s = merged["infected_death_sd"].to_numpy()

    grn_m = merged["altruistic_death_mean"].to_numpy()
    grn_s = merged["altruistic_death_sd"].to_numpy()

    gry_m = merged["cumulative_infection_mean"].to_numpy()
    gry_s = merged["cumulative_infection_sd"].to_numpy()

    fig, ax = plt.subplots(figsize=(4.2, 3.2), dpi=300)

    # ---- Cumulative infection (black)----
    lo, hi = band(gry_m, gry_s, n)
    ax.fill_between(T, lo, hi, color=COL_CUMUL, alpha=SHADE_ALPHA, linewidth=0)
    ax.plot(T, gry_m, color=COL_CUMUL, linewidth=2)

    # ---- Infected death (red) ----
    lo, hi = band(red_m, red_s, n)
    ax.fill_between(T, lo, hi, color=COL_INFECTED, alpha=SHADE_ALPHA, linewidth=0)
    ax.plot(T, red_m, color=COL_INFECTED, linewidth=2)

    # ---- Altruistic death (green) ----
    lo, hi = band(grn_m, grn_s, n)
    ax.fill_between(T, lo, hi, color=COL_ALTRUIST, alpha=SHADE_ALPHA, linewidth=0)
    ax.plot(T, grn_m, color=COL_ALTRUIST, linewidth=2)

    # ===== axis and format =====
    ax.set_xlabel("Time")
    ax.set_ylabel("population size")
    ax.set_xlim(T.min(), T.max())

    # No legend
    ax.legend_.remove() if ax.legend_ else None

    # remove some axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)

    print(f"Saved: {out_pdf}")
    
# ===== Execution =====
# run0_tstart0.csv ~ run39_tstart0.csv
df = load_runs("../../simulation/wt/out/run*_tstart_0.csv")
summ = summarize_counts(df)
merged = mean_and_band(summ)
plot_fig2c(merged, out_pdf="fig2c_repro.pdf", use_sem=False)
