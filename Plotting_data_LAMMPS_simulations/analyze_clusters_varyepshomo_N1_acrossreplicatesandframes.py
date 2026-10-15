import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from scipy.optimize import curve_fit
import argparse
import os
import sys
import csv
import seaborn as sns

sns.set()
sns.set_style("white")
sns.set_style("ticks")
sns.set_context("poster")

plt.rcParams.update({
    'font.family': "sans-serif",
    'font.size': 14,
    'axes.labelsize': 22,
    'axes.titlesize': 22,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'legend.frameon': False,
    'legend.framealpha': 0.8,
    'axes.linewidth': 1,
    'lines.linewidth': 2,
    'xtick.major.size': 8,
    'ytick.major.size': 8,
    'xtick.minor.size': 4,
    'ytick.minor.size': 4,
})

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def sigmoid(x, bottom, top, x0, k):
    return bottom + (top - bottom) / (1.0 + np.exp(-k * (x - x0)))


def compute_largest_cluster_fraction(positions, box, cutoff=2.0):
    n = len(positions)
    if n == 0:
        return 0.0, 0, 0

    box_lengths = box[:, 1] - box[:, 0]

    pos_wrapped = positions.copy()
    for dim in range(2):
        pos_wrapped[:, dim] = (
            (pos_wrapped[:, dim] - box[dim, 0]) % box_lengths[dim] + box[dim, 0]
        )

    images = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            shifted = pos_wrapped.copy()
            shifted[:, 0] += dx * box_lengths[0]
            shifted[:, 1] += dy * box_lengths[1]
            images.append(shifted)

    all_images = np.vstack(images)
    tree_orig = cKDTree(pos_wrapped)
    tree_images = cKDTree(all_images)
    pairs = tree_orig.query_ball_tree(tree_images, cutoff)

    rows, cols = [], []
    for i, neighbors in enumerate(pairs):
        for j_img in neighbors:
            j = j_img % n
            if i != j:
                rows.append(i)
                cols.append(j)

    if len(rows) == 0:
        return 1.0 / n, n, 1

    adj = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    n_components, labels = connected_components(adj, directed=False)
    cluster_sizes = np.bincount(labels)
    largest = np.max(cluster_sizes)

    return largest / n, n_components, largest


def load_frames(search_dir, args):
    if args.format == "lammpstrj":
        from lammps_traj_reader import read_lammpstrj
        all_frames = []
        for rep in range(1, args.n_reps + 1):
            rep_dir = os.path.join(search_dir, f"rep{rep}")
            trajpath = os.path.join(rep_dir, args.trajfile)
            if not os.path.exists(trajpath):
                print(f"WARNING: {trajpath} not found, skipping.")
                continue
            frames = read_lammpstrj(trajpath, atom_type=1,
                                    last_n_frames=args.last_n_frames,
                                    verbose=False)
            all_frames.extend(frames)
        return all_frames

    elif args.format == "datafile":
        from lammps_data_reader import read_lammps_data
        all_frames = []
        for rep in range(1, args.n_reps + 1):
            rep_dir = os.path.join(search_dir, f"rep{rep}")
            datapath = os.path.join(rep_dir, args.datafile)
            if not os.path.exists(datapath):
                continue
            frame = read_lammps_data(datapath, atom_type=1)
            if frame is not None:
                all_frames.append(frame)
        return all_frames

    return []


def main():
    parser = argparse.ArgumentParser(
        description="Cluster fraction vs eps_homo at fixed N")
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--values", nargs="+", type=float, required=True,
                        help="eps_homo values to sweep")
    parser.add_argument("--N_folder", default="N_1",
                        help="Fixed N folder (default: N_1)")
    parser.add_argument("--format", choices=["datafile", "lammpstrj"],
                        default="lammpstrj")
    parser.add_argument("--datafile", default="data_final_frame.xyz")
    parser.add_argument("--trajfile", default="file_trajectory.lammpstrj")
    parser.add_argument("--n_reps", type=int, default=5)
    parser.add_argument("--last_n_frames", type=int, default=10)
    parser.add_argument("--cutoff", type=float, default=2.0)
    parser.add_argument("--eps_hetero", type=float, default=15.0)
    parser.add_argument("--output", default="largestcluster_vs_epshomo")
    args = parser.parse_args()

    param_vals = []
    frac_means = []
    frac_stds = []

    for val in sorted(args.values):
        search_dir = os.path.join(args.base_dir, f"eps{val}", args.N_folder)

        if not os.path.isdir(search_dir):
            print(f"WARNING: {search_dir} not found, skipping.")
            continue

        print(f"eps_homo={val}:")
        frames = load_frames(search_dir, args)

        if len(frames) == 0:
            print(f"No frames found, skipping.")
            continue

        fracs = []
        for frame in frames:
            frac, nc, largest = compute_largest_cluster_fraction(
                frame["positions"], frame["box"], cutoff=args.cutoff
            )
            fracs.append(frac)

        param_vals.append(val)
        frac_means.append(np.mean(fracs))
        frac_stds.append(np.std(fracs))

        print(f"f_lc = {np.mean(fracs):.3f} +/- {np.std(fracs):.3f} "
              f"({len(frames)} frames)")

    if len(param_vals) == 0:
        print("No data found.")
        return

    param_vals = np.array(param_vals)
    frac_means = np.array(frac_means)
    frac_stds = np.array(frac_stds)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.errorbar(param_vals, frac_means, yerr=frac_stds, fmt="o",
                capsize=4, color="crimson", markersize=8, linewidth=1.5, zorder=5)

    try:
        p0 = [np.min(frac_means), np.max(frac_means),
              np.mean(param_vals), 10.0]
        bounds = ([0, 0, param_vals.min(), 0.1],
                  [1, 1, param_vals.max(), 200.0])
        popt, pcov = curve_fit(sigmoid, param_vals, frac_means,
                               p0=p0, bounds=bounds, maxfev=10000)
        x_fit = np.linspace(param_vals.min(), param_vals.max(), 200)
        y_fit = sigmoid(x_fit, *popt)
        ax.plot(x_fit, y_fit, "-", color="k", linewidth=1.3,
                label=f"Sigmoid fit", zorder=4)
        ax.axvline(x=popt[2], color="k", linestyle="--",
                   alpha=0.9, linewidth=1.2)
        print(f"\nSigmoid fit:")
        print(f"  bottom = {popt[0]:.3f}")
        print(f"  top    = {popt[1]:.3f}")
        print(f"  x0     = {popt[2]:.4f}  --> critical eps_homo")
        print(f"  k      = {popt[3]:.2f}")
    except RuntimeError:
        print("\n  WARNING: Sigmoid fit failed.")

    ax.set_xlabel(r"$\varepsilon_\mathrm{homo}$")
    ax.set_ylabel(r"$f_\mathrm{lc}$")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="best", fontsize=12)
    ax.set_title(rf"$\varepsilon_{{\mathrm{{hetero}}}}$ = {args.eps_hetero}")

    fig.tight_layout()
    fig.savefig(args.output + f"_epshetero{args.eps_hetero}.pdf", dpi=300)

if __name__ == "__main__":
    main()
