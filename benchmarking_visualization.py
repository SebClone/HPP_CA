# %%
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
os.makedirs("results/bench_results", exist_ok=True)
# %%
def plot_best_vs_baseline(data):
    """Plot best configuration runtime vs baseline for encrypt and decrypt modes in a single figure."""
    fig, ax = plt.subplots()
    for mode_name in ["encrypt", "decrypt"]:
        data_mode = data[data["mode"] == mode_name].copy()
        # Best per grid
        best = data_mode.loc[data_mode.groupby("grid")["runtime_s"].idxmin()].sort_values("grid")

        # Baseline per grid (np=1, omp_threads=1)
        baseline = data_mode[(data_mode["np"] == 1) & (data_mode["omp_threads"] == 1)].sort_values("grid")

        # Merge to ensure same grid ordering
        merged = pd.merge(
            best[["grid", "runtime_s"]],
            baseline[["grid", "runtime_s"]],
            on="grid",
            suffixes=("_best", "_baseline"),
        )

        if mode_name == "encrypt":
            ax.plot(merged["grid"], merged["runtime_s_best"], marker="o", label="Encrypt best", linestyle="-")
            ax.plot(merged["grid"], merged["runtime_s_baseline"], marker="o", linestyle="--", label="Encrypt baseline")
        else:
            ax.plot(merged["grid"], merged["runtime_s_best"], marker="^", label="Decrypt best", linestyle="-")
            ax.plot(merged["grid"], merged["runtime_s_baseline"], marker="^", linestyle="--", label="Decrypt baseline")

    ax.set_xlabel("Grid size")
    ax.set_ylabel("Runtime [s]")
    ax.set_title("Best config vs. baseline — encrypt & decrypt")
    ax.legend()
    ax.grid(True, linestyle=":")
    plt.tight_layout()
    fig.savefig("bench_results/best_vs_baseline_combined.png", dpi=300)


def strong_scaling(data, mode_name, grid_size, num_itrs):
    """Plot strong scaling results for a given mode and grid size."""
    data_str_scaling = data[(data["mode"] == mode_name) & (data["grid"] == grid_size)].copy()
    data_str_scaling["p"] = data_str_scaling["np"] * data_str_scaling["omp_threads"]

    # Baseline T1 (np=1, omp=1)
    base = data_str_scaling[(data_str_scaling["np"] == 1) & (data_str_scaling["omp_threads"] == 1)]["runtime_s"]
    if base.empty:
        print(f"[strong] missing baseline p=1 for mode={mode_name}, grid={grid_size}")
        return
    T1 = float(base.min())

    # Best runtime pro p
    best = data_str_scaling.groupby("p", as_index=False)["runtime_s"].min().sort_values("p")
    best["speedup"] = T1 / best["runtime_s"]

    # Plot
    fig, ax = plt.subplots()
    ax.plot(best["p"], best["speedup"], marker="o", label="Measured (best-of)")
    ax.plot(best["p"], best["p"].astype(float), linestyle=":", label="Ideal S=p")
    ax.set_xlabel("Total cores (np × omp_threads)")
    ax.set_ylabel("Speedup vs 1 core")
    ax.set_title(f"Strong scaling — {mode_name}, grid={grid_size}, iters={num_itrs}")
    ax.legend(); ax.grid(True, linestyle=":"); plt.tight_layout()
    out_path = f"bench_results/strong_simple_{mode_name}_g{grid_size}.png"
    fig.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")


def strong_scaling_combined(data, grids, num_itrs):
    """Plot strong scaling for multiple grids with encrypt (solid) and decrypt (dashed)."""
    fig, ax = plt.subplots()

    markers = ["o", "s", "D", "^", "v", "<", ">", "P", "X", "*"]
    marker_by_grid = {}

    # Determine a stable color per grid by plotting encrypt first to capture the color
    color_by_grid = {}

    # Build ideal p curve over the union of available p
    all_p = []
    for g in grids:
        df_tmp = data[(data["grid"] == g)].copy()
        if df_tmp.empty:
            continue
        df_tmp["p"] = df_tmp["np"] * df_tmp["omp_threads"]
        all_p.extend(df_tmp["p"].unique().tolist())
    if all_p:
        p_sorted = sorted(set(all_p))
        ax.plot(p_sorted, [float(p) for p in p_sorted], linestyle=":", label="Ideal S=p")

    for g in grids:
        # Ensure deterministic order: encrypt first, then decrypt to share color
        for mode_name in ["encrypt", "decrypt"]:
            df = data[(data["mode"] == mode_name) & (data["grid"] == g)].copy()
            if df.empty:
                print(f"[strong] no data for mode={mode_name}, grid={g}")
                continue

            df["p"] = df["np"] * df["omp_threads"]

            base = df[(df["np"] == 1) & (df["omp_threads"] == 1)]["runtime_s"]
            if base.empty:
                print(f"[strong] missing baseline p=1 for mode={mode_name}, grid={g}")
                continue
            T1 = float(base.min())

            best = df.groupby("p", as_index=False)["runtime_s"].min().sort_values("p")
            best["speedup"] = T1 / best["runtime_s"]

            linestyle = "-" if mode_name == "encrypt" else "--"

            if mode_name == "encrypt":
                marker = markers[len(marker_by_grid) % len(markers)]
                marker_by_grid[g] = marker
                line, = ax.plot(best["p"], best["speedup"], marker=marker, linestyle=linestyle, label=f"{mode_name} g={g}")
                color_by_grid[g] = line.get_color()
            else:
                marker = marker_by_grid.get(g, None)
                base_color = mcolors.to_rgb(color_by_grid[g])
                darker_color = tuple([c*0.7 for c in base_color])
                ax.plot(
                    best["p"], best["speedup"], marker=marker, markerfacecolor="none",
                    linestyle=linestyle, label=f"{mode_name} g={g}", color=darker_color
                )

    ax.set_xlabel("Total cores (np × omp_threads)")
    ax.set_ylabel("Speedup vs 1 core")
    ax.set_title(f"Strong scaling — encrypt (solid) + decrypt (dashed); iters={num_itrs}")
    ax.legend()
    ax.grid(True, linestyle=":")
    plt.tight_layout()
    out_path = "results/bench_results/strong_combined_allgrids.png"
    fig.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")


def weak_scaling(data, modes, base_grid, num_itrs):
    """Plot weak scaling for encrypt and decrypt in a single figure."""
    fig, ax = plt.subplots()
    for mode_name in modes:
        data_str_scaling = data[data["mode"] == mode_name].copy()
        data_str_scaling["p"] = data_str_scaling["np"] * data_str_scaling["omp_threads"]
        data_str_scaling["tau"] = (data_str_scaling["runtime_s"] / float(num_itrs)) * (data_str_scaling["p"] * (base_grid ** 2)) / (data_str_scaling["grid"].astype(float) ** 2)
        tau_best = data_str_scaling.groupby("p", as_index=False)["tau"].median().sort_values("p")
        linestyle = "-" if mode_name == "encrypt" else "--"
        if mode_name == "decrypt":
            ax.plot(tau_best["p"], tau_best["tau"], marker="o", markerfacecolor="none", linestyle=linestyle, label=mode_name)
        else:
            ax.plot(tau_best["p"], tau_best["tau"], marker="o", linestyle=linestyle, label=mode_name)
    ax.set_xlabel("Total cores (np × omp_threads)")
    ax.set_ylabel("Normalized time/iter (τ)")
    ax.set_title(f"Weak scaling (normalized), base_grid={base_grid}, iters={num_itrs}")
    ax.legend()
    ax.grid(True, linestyle=":")
    plt.tight_layout()
    out_path = f"results/bench_results/weak_simple_combined_base{base_grid}.png"
    fig.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")


# %%
num_itrs = 1000  # Fixed iterations for benchmarking
benchmark_data = pd.read_csv('results/bench_results/hpp_bench.csv')
benchmark_data = benchmark_data[benchmark_data['iters'] == num_itrs]
benchmark_data.drop(columns=['repeat', 'iters'], inplace=True)


encrypt_data = benchmark_data[benchmark_data['mode'] == 'encrypt']
decrypt_data = benchmark_data[benchmark_data['mode'] == 'decrypt']

# %%
# Best configuration per grid size
# Overall best per grid (independent of mode)
best_overall = benchmark_data.loc[
    benchmark_data.groupby('grid')['runtime_s'].idxmin()
].sort_values('grid').reset_index(drop=True)

# Best per (mode, grid)
best_by_mode = benchmark_data.loc[
    benchmark_data.groupby(['mode', 'grid'])['runtime_s'].idxmin()
].sort_values(['mode', 'grid']).reset_index(drop=True)

print("\nBest configuration per grid (overall):")
print(best_overall.to_string(index=False))

print("\nBest configuration per grid and mode:")
print(best_by_mode.to_string(index=False))

# Save to csv
best_overall.to_csv('results/bench_results/best_configs_overall.csv', index=False)
best_by_mode.to_csv('results/bench_results/best_configs_by_mode.csv', index=False)


# %%
# Plot best vs. baseline (1 process, 1 thread) for both modes in one figure
plot_best_vs_baseline(benchmark_data)

# %%
# Strong scaling kombiniert in einem Plot (Farben pro Grid, encrypt=solid, decrypt=dashed)
all_grids = sorted(benchmark_data["grid"].unique())
strong_scaling_combined(benchmark_data, all_grids, num_itrs)

# Weak scaling: nimm die kleinste Gridgröße als Referenz
base_grid_choice = int(sorted(benchmark_data["grid"].unique())[0])
weak_scaling(benchmark_data, ["encrypt", "decrypt"], base_grid_choice, num_itrs)

# %%
