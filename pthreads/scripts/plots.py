import argparse
from pathlib import Path
import re
from statistics import geometric_mean, stdev
from sys import implementation
from typing import List, Tuple
import sbatchman as sbm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools

### === Constants and Plot Setup === ###
FONT_TITLE = 30
FONT_AXES = 28
FONT_TICKS = 20
FONT_LEGEND = 18
COLOR_CYCLE = plt.rcParams["axes.prop_cycle"].by_key()["color"]
SET_FIG_TITLE = False
CHUNKSIZES = [4, 8, 16, 32, 64, 256, 1024, 4096]
OUT_DIR = Path("results")
OUT_DIR.mkdir(parents=True, exist_ok=True)
plt.style.use("seaborn-v0_8-colorblind")

plt.rc("axes", titlesize=FONT_AXES)
plt.rc("axes", labelsize=FONT_AXES)
plt.rc("xtick", labelsize=FONT_TICKS)
plt.rc("ytick", labelsize=FONT_TICKS)
plt.rc("legend", fontsize=FONT_LEGEND)
plt.rc("figure", titlesize=FONT_TITLE)

USE_LOG_SCALE = True

VALID_IMPLEMENTATIONS = ['atomic', 'mutex', 'openmp']
VALID_CHUNKSIZES = [4, 32, 256, 4096] # CHUNKSIZES
BEST_CHUNKSIZE = 32
IMPLEMENTATION_NAMES_MAP = {
  'atomic': 'Atomic',
  'mutex':  'Mutex',
  'openmp': 'OpenMP',
}
BOARD_NAMES_MAP = {
  'brah': 'AMD EPYC 7742',
  'baldo': 'AMD EPYC 7742',
  'pioneer': 'Milk-V Pioneer',
  'bananaf3': 'Banana Pi F3',
  'arriesgado': 'HiFive Unmatched',
}
BOARD_SHORT_NAMES_MAP = {
  'brah': 'AMD',
  'baldo': 'AMD',
  'AMD EPYC 7742': 'AMD',
  
  'pioneer': 'Pioneer',
  'Milk-V Pioneer': 'Pioneer',
  
  'bananaf3': 'BananaPi',
  'Banana Pi F3': 'BananaPi',
  
  'arriesgado': 'HiFive',
  'HiFive Unmatched': 'HiFive',
}

### === Functions === ###
def filter_jobs(jobs: List[sbm.Job]) -> List[sbm.Job]:
  return [job for job in jobs if job.status == "COMPLETED" and "compile" not in job.tag]


def parse_stdout(job: sbm.Job) -> Tuple[float, float, float, float]:
  stdout = job.get_stdout()
  times = []
  for line in stdout.splitlines():
    if line.startswith("run_id=") or line.startswith("run="):
      try:
        times.append(float(line.split(",")[-1]))
      except ValueError:
        continue
  times = times[2:] # Warmup
  return geometric_mean(times) if times else np.nan, stdev(times), min(times), max(times)


def jobs_to_dataframe() -> pd.DataFrame:
  jobs = sbm.jobs_list(from_active=True, from_archived=True)
  jobs = filter_jobs(jobs)

  parsed_jobs = []
  for job in jobs:
    _, __, args = job.parse_command_args()
    avg_runtime, std_runtime, min_runtime, max_runtime = parse_stdout(job)

    if job.tag.startswith("openmp"):
      chunksizes = CHUNKSIZES
      implementation = "openmp"
      dataset = Path(job.command.split()[1]).stem
    else:
      chunksizes = [int(re.match(r"\w+chunksize_(\d+)", job.tag).group(1))]
      implementation = re.match(r"(\w+)chunksize_\d+", job.tag).group(1)[:-1]
      dataset = Path(args["f"]).stem

    for chunksize in chunksizes:
      parsed_jobs.append({
        "dataset": dataset,
        "board": sbm.get_cluster_name(),
        "num_cpus": int(re.match(r"\w*(\d+)_?cpus", job.config_name).group(1)),
        "chunksize": chunksize,
        "implementation": implementation,
        "runtime": avg_runtime,
        "std_runtime": std_runtime,
        "min_runtime": min_runtime,
        "max_runtime": max_runtime,
      })

  df = pd.DataFrame(parsed_jobs)
  df.dropna(inplace=True)
  out_path = OUT_DIR / "data.csv"
  df.to_csv(out_path, index=False)
  print(f"CSV data saved to {out_path.resolve()}")
  return df


def load_csv_files(filepaths: List[str]) -> pd.DataFrame:
  dfs = []
  for file in filepaths:
    cluster = Path(file).stem.split("_")[0]
    df = pd.read_csv(file)
    df["cluster"] = cluster
    dfs.append(df)
  return pd.concat(dfs, ignore_index=True)


def line_plot(ax, data, title, color_map, marker='o'):
  for impl, grp in data.groupby("implementation"):
    mean_grp = grp.groupby("num_cpus", as_index=False)["runtime"].mean()
    ax.plot(
      mean_grp["num_cpus"],
      mean_grp["runtime"],
      marker=marker,
      label=IMPLEMENTATION_NAMES_MAP.get(impl, impl),
      color=color_map.get(impl),
    )
  if USE_LOG_SCALE: ax.set_xscale('log', base=2)
  cores_unique = sorted(data['num_cpus'].unique())
  ax.set_xticks(cores_unique)
  ax.set_xticklabels(cores_unique)
  ax.set_xlabel("Cores")
  ax.set_ylabel("Runtime [s]")
  ax.set_title(title, fontsize=FONT_TITLE-10)
  ax.grid(True)
  ax.legend() # title="Implementation")


def generate_line_plots(df: pd.DataFrame):
  implementations = df["implementation"].unique()
  color_map = dict(zip(implementations, itertools.cycle(COLOR_CYCLE)))

  for board_name, df_board in df.groupby("board"):
    for dataset_name, df_dataset in df_board.groupby("dataset"):
      y_min, y_max = df_dataset["runtime"].min() * 0.95, df_dataset["runtime"].max() * 1.05
      df_dataset = df_dataset.sort_values("num_cpus")
      x_ticks = sorted(df_dataset["num_cpus"].unique())
      chunks = sorted(df_dataset["chunksize"].unique())
      n_cols = int(len(chunks)**0.5 + 0.999)
      n_rows = (len(chunks) + n_cols - 1) // n_cols

      fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
      axes = axes.flat if len(chunks) > 1 else [axes]

      for ax, ch in zip(axes, chunks):
        ax.set_ylim(y_min, y_max)
        ax.set_xticks(x_ticks)
        line_plot(ax, df_dataset[df_dataset["chunksize"] == ch], f"Chunksize = {ch}", color_map)

      for ax in axes[len(chunks):]:
        ax.axis("off")

      if SET_FIG_TITLE:
        fig.suptitle(f"{board_name} - Strong Scaling (graph: {dataset_name})", y=0.98)
        
      fig.tight_layout()
      path = OUT_DIR / 'scaling' / f'{board_name.replace(" ", "_")}_strong_scaling_{dataset_name.replace(" ", "_")}.png'
      path.parent.mkdir(exist_ok=True, parents=True)
      plt.savefig(path)
      plt.close()
      print(f"Plot saved to {path.resolve()}")


def generate_line_plots_for_board_pairs(df: pd.DataFrame, board_pairs: list[tuple[str, str]]):
  """
  Generate combined strong scaling plots for each dataset, showing data from both boards in the given pairs.
  Each subplot corresponds to a chunksize; lines are colored by implementation and styled by board.
  """
  implementations = sorted(df["implementation"].unique())
  color_map = dict(zip(sorted(implementations), itertools.cycle(COLOR_CYCLE)))

  # Distinct linestyles per board
  line_styles = ['-', '--', '-.', ':']
  board_linestyle_map = {}
  all_boards = set(itertools.chain.from_iterable(board_pairs))
  for style, board in zip(itertools.cycle(line_styles), sorted(all_boards)):
    board_linestyle_map[board] = style

  for board1, board2 in board_pairs:
    df_pair = df[df["board"].isin([board1, board2])]
    for dataset_name, df_dataset in df_pair.groupby("dataset"):
      y_min, y_max = df_dataset["runtime"].min() * 0.95, df_dataset["runtime"].max() * 1.05
      chunks = sorted(df_dataset["chunksize"].unique())
      n_cols = int(len(chunks) ** 0.5 + 0.999)
      n_rows = (len(chunks) + n_cols - 1) // n_cols

      fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4.5 * n_rows))
      axes = axes.flat if len(chunks) > 1 else [axes]

      for idx, (ax, ch) in enumerate(zip(axes, chunks)):
        # UNCOMMENT TO GET SAME SCALE ACROSS PLOTS ax.set_ylim(y_min, y_max)
        ax.set_title(f"Chunksize = {ch}", fontsize=FONT_TITLE-10)
        if USE_LOG_SCALE: ax.set_xscale('log', base=2)

        for (board, impl), group in df_dataset[(df_dataset["chunksize"] == ch)].groupby(["board", "implementation"]):
          label = f"{IMPLEMENTATION_NAMES_MAP.get(impl, impl)} ({BOARD_SHORT_NAMES_MAP.get(board,board)})"
          linestyle = board_linestyle_map[board]
          color = color_map[impl]
          group = group.sort_values(["num_cpus"])
          mean_grp = group.groupby("num_cpus", as_index=False)["runtime"].mean()
          ax.plot(mean_grp["num_cpus"], mean_grp["runtime"], label=label, color=color, linestyle=linestyle, marker='o', markersize=3)

        cores_unique = sorted(df_dataset['num_cpus'].unique())
        ax.set_xticks(cores_unique)
        ax.set_xticklabels(cores_unique)
        ax.grid(True)
        ax.legend(frameon=False)

        if idx // n_cols == n_rows - 1:
          ax.set_xlabel("Cores")
        if idx % n_cols == 0:
          ax.set_ylabel("Runtime [s]")

      handles, labels = axes[0].get_legend_handles_labels()
      sorted_items = sorted(zip(labels, handles), key= lambda x: x[0].split(" ")[0]) # x[0].split("(")[-1].strip(")"))
      labels, handles = zip(*sorted_items) if sorted_items else ([], [])
      fig.legend(
        handles, labels,
        loc='upper center',
        ncol=len(labels)/2,
        frameon=False
      )

      for ax in axes:
        # ax.axis("off")
        if ax.get_legend() is not None:
          ax.get_legend().remove()
        
      if SET_FIG_TITLE:
        fig.suptitle(f"{board1} vs {board2} - Strong Scaling (graph: {dataset_name})", y=0.98)
      fig.tight_layout(rect=[0., 0., 1, 0.92])
      path = OUT_DIR / 'scaling_comparison' / f'{board1.replace(" ", "_")}_vs_{board2.replace(" ", "_")}_strong_scaling_{dataset_name.replace(" ", "_")}.png'
      path.parent.mkdir(exist_ok=True, parents=True)
      plt.savefig(path)
      plt.close()
      print(f"Plot saved to {path.resolve().absolute()}")


def plot_board_comparisons(df: pd.DataFrame, min_coverage_ratio=0.5):
  df = df[df['chunksize'] == BEST_CHUNKSIZE]
  boards = df["board"].unique()
  datasets = df["dataset"].unique()
  implementations = sorted(df["implementation"].unique())

  # Assign unique colors to implementations
  colors = list(COLOR_CYCLE)
  impl_colors = {impl: colors[i] for i, impl in enumerate(sorted(implementations))}

  for dataset in datasets:
    df_dataset = df[df["dataset"] == dataset]

    for board1, board2 in [(BOARD_NAMES_MAP['brah'], BOARD_NAMES_MAP['pioneer'])]: # FIXME itertools.combinations(boards, 2):
      fig, ax = plt.subplots(figsize=(12, 5))
      bar_width = 0.8 / len(implementations)
      miny, maxy = np.inf, -np.inf

      # Build merged data for all implementations
      cpu_data = {}
      for impl in implementations:
        impl_df = df_dataset[df_dataset["implementation"] == impl]
        b1_df = impl_df[impl_df["board"] == board1]
        b2_df = impl_df[impl_df["board"] == board2]
        merged = pd.merge(
          b1_df, b2_df,
          on="num_cpus",
          suffixes=(f"_{board1}", f"_{board2}")
        )
        if not merged.empty:
          r1 = merged[f"runtime_{board1}"]
          r2 = merged[f"runtime_{board2}"]
          # merged["norm_diff"] = (r1 - r2) / r1
          merged["speedup"] = np.where(r1 < r2, -(r2 / r1) + 1, (r1 / r2) - 1)
          
          for _, row in merged.iterrows():
            cpu = row["num_cpus"]
            if cpu not in cpu_data:
              cpu_data[cpu] = {}
            cpu_data[cpu][impl] = row["speedup"] # (row["norm_diff"], row["speedup"])
        # if impl == 'mutex':
        #   print(merged[['num_cpus', 'runtime_AMD EPYC 7742', 'runtime_Milk-V Pioneer','speedup']].sort_values('num_cpus'))

      # Filter CPU counts by coverage
      valid_cpu_counts = [
        cpu for cpu, impls in cpu_data.items()
        if len(impls) / len(implementations) >= min_coverage_ratio
      ]

      if not valid_cpu_counts:
        plt.close()
        continue

      valid_cpu_counts = sorted(valid_cpu_counts)
      x_base = np.arange(len(valid_cpu_counts))
      offsets = np.linspace(-0.4 + bar_width / 2, 0.4 - bar_width / 2, len(implementations))

      for i, impl in enumerate(implementations):
        # norm_diffs = []
        speedups = []
        x_pos = []

        for idx, cpu in enumerate(valid_cpu_counts):
          if impl in cpu_data[cpu]:
            spd = cpu_data[cpu][impl]
            # norm_diffs.append(norm)
            speedups.append(spd)
            x_pos.append(x_base[idx] + offsets[i])

        if not speedups:
          continue

        bars = ax.bar(
          x_pos,
          speedups,
          width=bar_width,
          color=impl_colors[impl],
          label=IMPLEMENTATION_NAMES_MAP.get(impl, impl)
        )
        ax.axhline(0, color='r', linewidth=3)

        for bar, speedup in zip(bars, speedups):
          height = bar.get_height()
          offset = 0.04
          label_y = height + offset if height >= 0 else height - offset
          va = "bottom" if height >= 0 else "top"
          ax.text(bar.get_x() + bar.get_width() / 2, label_y, f"{abs(speedup)+1:.2f}x", ha="center", va=va, fontsize=10)
          miny, maxy = min(miny, speedup), max(maxy, speedup)

      ax.axhline(0, color="black", linewidth=1)
      ax.set_xticks(x_base)
      ax.set_xticklabels([str(cpu) for cpu in valid_cpu_counts])
      yticks = np.arange(-np.ceil(np.abs(miny)), np.ceil(maxy), 1)
      ax.set_yticks(yticks)
      yticklabels = []
      yticklabels_start = yticks[0]-1
      for i in range(len(yticks)):
        tick = yticklabels_start+i
        if tick >= -1:
          tick += 2
        yticklabels.append(f'${int(abs(tick))}\\times$') # if int(tick)!= 1 else 'Equal')
      ax.set_yticklabels(yticklabels)
      ax.set_xlabel("Cores")
      ax.set_ylabel("Speedup")
      if SET_FIG_TITLE:
        ax.set_title(f"{board1} vs {board2} (Dataset: {dataset}, ChunkSize: {BEST_CHUNKSIZE})")
      miny = yticks[0]-1.0
      maxy = yticks[-1]+1.0
      ax.set_ylim(miny, maxy)
      ax.text(-0.7, maxy-0.1,  f'{board2} Faster', fontsize=18, horizontalalignment='left', verticalalignment='top')
      ax.text(-0.7, miny+0.05, f'{board1} Faster', fontsize=18, horizontalalignment='left', verticalalignment='bottom')
      ax.legend() # title="Implementation")
      ax.grid(True, axis="y", linestyle="--", alpha=0.6)
      fig.tight_layout()

      out_path = OUT_DIR / 'comparison' / f'{dataset.replace(" ", "_")}_compare_{board1.replace(" ", "_")}_vs_{board2.replace(" ", "_")}_grouped.png'
      out_path.parent.mkdir(exist_ok=True, parents=True)
      plt.savefig(out_path)
      plt.close()
      print(f"Grouped comparison plot saved to {out_path.resolve()}")


### === Main CLI === ###
def main():
  parser = argparse.ArgumentParser(description="Plot BFS benchmark results.")
  parser.add_argument("csv_files", nargs="*", help="CSV file(s) containing BFS results")
  args = parser.parse_args()

  if args.csv_files:
    df = load_csv_files(args.csv_files)
    print(f"Loaded data from {len(args.csv_files)} CSV file(s)")
  else:
    df = jobs_to_dataframe()

  # Filter
  df = df[df['implementation'].isin(VALID_IMPLEMENTATIONS)]
  df = df[df['chunksize'].isin(VALID_CHUNKSIZES)]
  # Remap board names
  df['board'] = df['board'].map(BOARD_NAMES_MAP)
  
  # generate_line_plots(df)
  # generate_line_plots_for_board_pairs(df, [(BOARD_NAMES_MAP["baldo"], BOARD_NAMES_MAP["pioneer"])])
  # plot_board_comparisons(df)


if __name__ == "__main__":
  main()
