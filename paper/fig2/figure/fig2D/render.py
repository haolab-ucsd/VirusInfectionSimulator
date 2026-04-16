import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_condition(pattern, tstart_label):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {pattern}")

    dfs = []
    for p in paths:
        m = re.search(r"run(\d+)_.*_tstart_(\d+)\.csv", p)
        run = int(m.group(1)) if m else None

        df = pd.read_csv(p, usecols=["T", "T_step", "S"])
        df["run"] = run
        df["tstart"] = tstart_label
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)

def summarize_population(df):
    g = df.groupby(["tstart", "run", "T_step", "T"])["S"]

    summ = g.agg(
        population=lambda s: np.sum(s.isin([1, 2, 4]))
    ).reset_index()

    return summ

def mean_over_runs(summ):
    grp = summ.groupby(["tstart", "T_step", "T"])

    mean = grp["population"].mean()
    sd   = grp["population"].std(ddof=1)
    n    = grp.size()

    out = pd.concat([mean, sd, n], axis=1).reset_index()
    out.columns = ["tstart", "T_step", "T", "mean", "sd", "n"]
    return out

def plot_figD(df_mean, out_pdf="figD_repro.pdf", use_sem=False):
    fig, ax = plt.subplots(figsize=(4.2, 3.2), dpi=300)

    colors = {
        "0":   "#000000",  # No pretreatment
        "50":  "#666666",
        "250": "#bbbbbb",
    }

    labels = {
        "0":   "No pretreatment",
        "50":  "t = 50 pretreatment",
        "250": "t = 250 pretreatment",
    }

    for tstart, sub in df_mean.groupby("tstart"):
        T = sub["T"].to_numpy()
        y = sub["mean"].to_numpy()

        if use_sem:
            band = sub["sd"].to_numpy() / np.sqrt(sub["n"].to_numpy())
        else:
            band = sub["sd"].to_numpy()

        ax.plot(
            T, y,
            color=colors[str(tstart)],
            linewidth=2,
            label=labels[str(tstart)]
        )

    ax.set_xlabel("Time Post Infection")
    ax.set_ylabel("Mean population size")
    ax.set_xlim(df_mean["T"].min(), df_mean["T"].max())

    # ax.legend(frameon=False, loc="upper left")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


df0   = load_condition("../../simulation/wt/out/run*_tstart_0.csv",   tstart_label="0")
df50  = load_condition("../../simulation/pretreatment_t50/out/run*_tstart_50.csv",  tstart_label="50")
df250 = load_condition("../../simulation/pretreatment_t250/out/run*_tstart_250.csv", tstart_label="250")

df_all = pd.concat([df0, df50, df250], ignore_index=True)


summ = summarize_population(df_all)
mean_df = mean_over_runs(summ)

plot_figD(mean_df, out_pdf="figD_repro.pdf", use_sem=False)
