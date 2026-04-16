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
        m = re.search(r"run(\d+)_tstart[_]?(\d+)\.csv", p)
        run = int(m.group(1)) if m else None

        df = pd.read_csv(p, usecols=["T", "T_step", "S"])
        df["run"] = run if run is not None else p
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def summarize_counts(df):
    """
    Updated version:
      - infected_death = DA + DE  (S in {2,4})
      - altruistic_death = DB     (S == 3)
      - cumulative_infection = I + DA + DE (S in {1,2,4})
    """
    g = df.groupby(["run", "T_step", "T"])["S"]

    out = g.agg(
        infected_death=lambda s: np.sum(s.isin([2, 4])),
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

    merged = mean.merge(sd, on=["T_step", "T"], suffixes=("_mean", "_sd")).merge(n, on=["T_step", "T"])
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
    COL_INFECTED = "#d62728"
    COL_ALTRUIST = "#2ca02c"
    COL_CUMUL    = "#000000"

    SHADE_ALPHA = 0.25

    red_m = merged["infected_death_mean"].to_numpy()
    red_s = merged["infected_death_sd"].to_numpy()

    grn_m = merged["altruistic_death_mean"].to_numpy()
    grn_s = merged["altruistic_death_sd"].to_numpy()

    gry_m = merged["cumulative_infection_mean"].to_numpy()
    gry_s = merged["cumulative_infection_sd"].to_numpy()

    fig, ax = plt.subplots(figsize=(4.2, 3.2), dpi=300)

    # ---- Cumulative infection (black) ----
    lo, hi = band(gry_m, gry_s, n)
    ax.fill_between(T, lo, hi, color=COL_CUMUL, alpha=SHADE_ALPHA, linewidth=0)
    ax.plot(T, gry_m, color=COL_CUMUL, linewidth=2)

    # ---- Infected death (red): now DA + DE ----
    lo, hi = band(red_m, red_s, n)
    ax.fill_between(T, lo, hi, color=COL_INFECTED, alpha=SHADE_ALPHA, linewidth=0)
    ax.plot(T, red_m, color=COL_INFECTED, linewidth=2)

    # ---- Altruistic death (green): DB ----
    lo, hi = band(grn_m, grn_s, n)
    ax.fill_between(T, lo, hi, color=COL_ALTRUIST, alpha=SHADE_ALPHA, linewidth=0)
    ax.plot(T, grn_m, color=COL_ALTRUIST, linewidth=2)

    ax.set_xlabel("Time")
    ax.set_ylabel("population size")
    ax.set_xlim(T.min(), T.max())

    ax.legend_.remove() if ax.legend_ else None
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)

    print(f"Saved: {out_pdf}")


def summarize_ratio_db_over_infected_death(df):
    """
    Updated version:
      - infected death denominator = DA + DE
      - altruistic / infected death = DB / (DA + DE)

    To mimic the earlier MATLAB-style handling of undefined ratios:
      - replace non-finite values with 0
    """
    g = df.groupby(["run", "T_step", "T"])["S"]

    out = g.agg(
        DA_count=lambda s: np.sum(s == 2),
        DB_count=lambda s: np.sum(s == 3),
        DE_count=lambda s: np.sum(s == 4),
        cumulative_infection=lambda s: np.sum(s.isin([1, 2, 4])),
    ).reset_index()

    denom = out["DA_count"] + out["DE_count"]

    with np.errstate(divide="ignore", invalid="ignore"):
        out["death_ratio"] = out["DB_count"] / denom

    out.loc[~np.isfinite(out["death_ratio"]), "death_ratio"] = 0.0

    return out


def mean_and_band_ratio(summ):
    grp = summ.groupby(["T_step", "T"])

    mean = grp[["cumulative_infection", "death_ratio"]].mean().reset_index()
    sd   = grp[["cumulative_infection", "death_ratio"]].std(ddof=1).reset_index()
    n    = grp.size().reset_index(name="n")

    merged = mean.merge(sd, on=["T_step", "T"], suffixes=("_mean", "_sd")).merge(n, on=["T_step", "T"])
    return merged


def plot_fig2c_like_fig1d(merged, out_pdf="fig2c_like_fig1d_DBoverDAplusDE.pdf", use_sem=False):
    T = merged["T"].to_numpy()
    n = merged["n"].to_numpy()

    def band(y_mean, y_sd, n):
        if use_sem:
            y_band = y_sd / np.sqrt(n)
        else:
            y_band = y_sd
        return y_mean - y_band, y_mean + y_band

    # ===== Color definition (Fig 1D-like) =====
    COL_CUMUL_LINE = "#3b4cc0"
    COL_CUMUL_FILL = "#9aa7ff"

    COL_RATIO_LINE = "#e07b39"
    COL_RATIO_FILL = "#f3b38a"

    SHADE_ALPHA = 0.35

    cumul_m = merged["cumulative_infection_mean"].to_numpy()
    cumul_s = merged["cumulative_infection_sd"].to_numpy()

    ratio_m = merged["death_ratio_mean"].to_numpy()
    ratio_s = merged["death_ratio_sd"].to_numpy()

    fig, ax1 = plt.subplots(figsize=(4.6, 3.4), dpi=300)
    ax2 = ax1.twinx()

    # ---- Left axis: cumulative infection (blue) ----
    lo, hi = band(cumul_m, cumul_s, n)
    ax1.fill_between(T, lo, hi, color=COL_CUMUL_FILL, alpha=SHADE_ALPHA, linewidth=0)
    ax1.plot(T, cumul_m, color=COL_CUMUL_LINE, linewidth=2)

    # ---- Right axis: DB / (DA + DE) (orange) ----
    lo, hi = band(ratio_m, ratio_s, n)
    ax2.fill_between(T, lo, hi, color=COL_RATIO_FILL, alpha=SHADE_ALPHA, linewidth=0)
    ax2.plot(T, ratio_m, color=COL_RATIO_LINE, linewidth=2)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Cumulative infection", color=COL_CUMUL_LINE)
    ax2.set_ylabel("Altruistic / Infection death", color=COL_RATIO_LINE)

    ax1.set_xlim(T.min(), T.max())

    ax1.tick_params(axis="y", colors=COL_CUMUL_LINE)
    ax2.tick_params(axis="y", colors=COL_RATIO_LINE)

    ax1.spines["left"].set_color(COL_CUMUL_LINE)
    ax2.spines["right"].set_color(COL_RATIO_LINE)

    ax1.yaxis.label.set_color(COL_CUMUL_LINE)
    ax2.yaxis.label.set_color(COL_RATIO_LINE)

    ymax = np.nanmax(ratio_m + ratio_s)
    if not np.isfinite(ymax) or ymax <= 0:
        ymax = 1.0
    ax2.set_ylim(0, ymax * 1.1)

    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)

    print(f"Saved: {out_pdf}")


# ===== Execution =====
df = load_runs("../../simulation/wt/out/run*_tstart_0.csv")

# Updated Fig 2C-like plot
summ = summarize_counts(df)
merged = mean_and_band(summ)
plot_fig2c(merged, out_pdf="fig2c_repro_DAplusDE.pdf", use_sem=False)

# Updated Fig 1D-like plot
summ_ratio = summarize_ratio_db_over_infected_death(df)
merged_ratio = mean_and_band_ratio(summ_ratio)
plot_fig2c_like_fig1d(
    merged_ratio,
    out_pdf="fig2c_like_fig1d_DBoverDAplusDE.pdf",
    use_sem=False
)
