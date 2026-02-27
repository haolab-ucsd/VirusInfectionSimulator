import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.collections import PatchCollection

# =========================
# Setting of the colour
# =========================
COL_UNINFECTED = "#c7b6ff"
COL_INFECTED   = "#ff2d2d"
COL_DEAD_U     = "#2fbf71"
COL_OTHER      = "#ffd84d"

# S: 0=UI, 1=I, 2=DA, 3=DB, 4=DE
STATE_TO_COLOR = {
    0: COL_UNINFECTED,
    1: COL_INFECTED,
    2: COL_OTHER,
    3: COL_DEAD_U,
    4: COL_OTHER,
}

def estimate_hex_radius(H: np.ndarray, V: np.ndarray, seed: int = 0) -> float:
    rng = np.random.default_rng(seed)
    pts = np.column_stack([H, V])
    n = len(pts)
    m = min(n, 5000)
    idx = rng.choice(n, size=m, replace=False)
    samp = pts[idx]

    dmin = []
    dx_all = np.abs(pts[:, 0] - np.median(pts[:, 0]))
    x_band = np.quantile(dx_all, 0.01) if n > 10 else 1.0

    for p in samp:
        dx = np.abs(pts[:, 0] - p[0])
        cand = pts[dx < x_band]
        if len(cand) < 2:
            continue
        d = np.sqrt(((cand - p) ** 2).sum(axis=1))
        d = d[d > 0]
        if len(d) > 0:
            dmin.append(d.min())

    if len(dmin) == 0:
        uniq = np.unique(np.round(H, 8))
        if len(uniq) >= 2:
            d_nn = np.median(np.diff(np.sort(uniq)))
        else:
            d_nn = 1.0
    else:
        d_nn = float(np.median(dmin))

    r_pointy = d_nn / np.sqrt(3)
    r_flat   = 2 * d_nn / 3
    r = 0.99 * min(r_pointy, r_flat)
    return float(r)

def make_all_frames_from_csv(
    csv_path: str,
    out_dir: str = "frames",
    time_col: str = "T_step",
    every: int = 1,
    dpi: int = 150,
    figsize: tuple = (7, 7),
    show_title: bool = True,
    orientation: float = 0.0,
    pad_factor: float = 2.2,
):
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(csv_path)

    required = {"H", "V", "S", time_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"No required data in input CSV: {missing}")

    times = np.sort(df[time_col].unique())
    times = times[::every]

    t0 = times[0]
    df0 = df[df[time_col] == t0].copy()

    H0 = df0["H"].to_numpy()
    V0 = df0["V"].to_numpy()

    r = estimate_hex_radius(H0, V0, seed=0)

    patches = [
        RegularPolygon((h, v), numVertices=6, radius=r, orientation=orientation)
        for h, v in zip(H0, V0)
    ]

    pad = pad_factor * r
    x_min, x_max = H0.min() - pad, H0.max() + pad
    y_min, y_max = V0.min() - pad, V0.max() + pad

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    S0 = df0["S"].to_numpy().astype(int)
    facecolors = [STATE_TO_COLOR.get(int(s), COL_OTHER) for s in S0]

    coll = PatchCollection(
        patches,
        facecolor=facecolors,
        edgecolor="none",
        linewidth=0,
        antialiased=False,
    )
    ax.add_collection(coll)

    for i, t in enumerate(times):
        sub = df[df[time_col] == t]

        if "IDX" in sub.columns and "IDX" in df0.columns:
            sub = sub.sort_values("IDX")
            if i == 0:
                df0_sorted = df0.sort_values("IDX")

        S = sub["S"].to_numpy().astype(int)
        colors = [STATE_TO_COLOR.get(int(s), COL_OTHER) for s in S]
        coll.set_facecolor(colors)

        if show_title:
            if "T" in sub.columns:
                Tval = float(sub["T"].iloc[0])
                ax.set_title(f"T = {Tval:.2f} ({time_col}={int(t)})", fontsize=12)
            else:
                ax.set_title(f"{time_col}={int(t)}", fontsize=12)

        out_path = os.path.join(out_dir, f"frame_{i:04d}.png")
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0.15)

        if (i + 1) % 10 == 0:
            print(f"[{i+1}/{len(times)}] saved: {out_path}")

    plt.close(fig)
    print(f"Done. Saved {len(times)} frames to '{out_dir}'")
    print("Next (Making movie!): ffmpeg -framerate 30 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p out.mp4")

if __name__ == "__main__":
    # ex :
    # python make_frames.py run0_tstart0.csv frames
    import sys
    csv_path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "frames"
    make_all_frames_from_csv(csv_path, out_dir=out_dir, every=1, dpi=150)
