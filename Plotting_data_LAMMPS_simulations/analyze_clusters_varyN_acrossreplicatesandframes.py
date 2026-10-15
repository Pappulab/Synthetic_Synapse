import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from scipy.optimize import curve_fit
import argparse
import os
import sys
import glob
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


def find_datafile(directory, pattern="data*.xyz"):
    matches = glob.glob(os.path.join(directory, pattern))
    if len(matches) == 0:
        return None
    return max(matches, key=os.path.getmtime)


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

    rows = []
    cols = []
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


def analyze_single_eps(base_dir, N_values, args):
    results = {
        "N": [],
        "cluster_frac_mean": [],
        "cluster_frac_std": []}

    for N in sorted(N_values):
        search_dir = os.path.join(base_dir, f"N_{N}")

        if not os.path.isdir(search_dir):
            print(f"WARNING: {search_dir} not found, skipping.")
            continue

        frames = load_frames(search_dir, args)
        if len(frames) == 0:
            print(f"WARNING: No frames in {search_dir}, skipping.")
            continue

        fracs = []

        for frame in frames:
            pos = frame["positions"]
            box = frame["box"]

            frac, nc, largest = compute_largest_cluster_fraction(
                pos, box, cutoff=args.cutoff
            )

            fracs.append(frac)

        results["N"].append(N)
        results["cluster_frac_mean"].append(np.mean(fracs))
        results["cluster_frac_std"].append(np.std(fracs))

        print(f"N={N:6d} | cluster: {np.mean(fracs):.3f} "
              f"+/- {np.std(fracs):.3f} | "
              f"({len(frames)} frames)")

    for key in results:
        results[key] = np.array(results[key])

    return results

def main():
    parser = argparse.ArgumentParser(
        description="Cluster analysis: N sweep for different eps_homo values")
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--eps_values", nargs="+", type=float, default=None,
                        help="Multiple eps_homo values to overlay")
    parser.add_argument("--values", nargs="+", type=int, required=True,
                        help="N (receptor count) values to sweep")
    parser.add_argument("--format", choices=["binary", "datafile", "lammpstrj"],
                        default="lammpstrj")
    parser.add_argument("--trajfile", default="file_trajectory.lammpstrj")
    parser.add_argument("--n_reps", type=int, default=5)
    parser.add_argument("--last_n_frames", type=int, default=10)
    parser.add_argument("--cutoff", type=float, default=2.0)
    parser.add_argument("--eps_hetero", type=float, default=50.0)
    parser.add_argument("--fit", action="store_true", default=False)
    parser.add_argument("--logx", action="store_true", default=False)
    parser.add_argument("--output", default="largestclusterscaffold_vs_N")
    args = parser.parse_args()

    N_values = sorted(args.values)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    if args.eps_values and len(args.eps_values) > 1:
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(args.eps_values)))
    else:
        colors = ["crimson"]

    all_results = {}  

    if args.eps_values:
        for i, eps in enumerate(sorted(args.eps_values)):
            eps_dir = os.path.join(args.base_dir, f"eps{eps}")
            if not os.path.isdir(eps_dir):
                print(f"  WARNING: {eps_dir} not found, skipping.")
                continue

            print(f"  eps_homo = {eps}:")
            results = analyze_single_eps(eps_dir, N_values, args)

            if len(results["N"]) == 0:
                continue

            all_results[eps] = results

            label = rf"$\varepsilon_{{\mathrm{{homo}}}}$ = {eps}"

            ax.errorbar(results["N"], results["cluster_frac_mean"],
                        yerr=results["cluster_frac_std"],
                        fmt="o-", color=colors[i], markersize=8,
                        linewidth=2, capsize=3, label=label)

    else:
        print(f"  Analyzing {args.base_dir}:")
        results = analyze_single_eps(args.base_dir, N_values, args)

        if len(results["N"]) == 0:
            print("No data found. Check paths.")
            return

        all_results["single"] = results  

        ax.errorbar(results["N"], results["cluster_frac_mean"],
                    yerr=results["cluster_frac_std"],
                    fmt="o", color="crimson", markersize=8, capsize=3)

    ax.set_xlabel(r"$N_\mathrm{receptor}$")
    ax.set_ylabel(r"$f_\mathrm{lc}$")
    ax.set_ylim(-0.05, 1.05)
    if args.logx:
        ax.set_xscale("log")
    if args.eps_values and len(args.eps_values) > 1:
        ax.legend(loc="best", fontsize=12)
    ax.set_title(rf"$\varepsilon_{{\mathrm{{hetero}}}}$ = {args.eps_hetero}")

    fig.tight_layout()
    fig.savefig(args.output + f"_epshetero{args.eps_hetero}.pdf", dpi=300)

if __name__ == "__main__":
    main()




