import glob
import os
from typing import Optional, Tuple
import pandas as pd
import numpy as np
import re
from statistics import stdev, geometric_mean

def parse_log_file(log_file_path: str) -> Tuple[float, float, float, float]:
  """Parses a log file and returns performance metrics."""
  try:
    with open(log_file_path, "r") as f:
      stdout = f.read()
  except FileNotFoundError:
    return np.nan, np.nan, np.nan, np.nan

  if not stdout:
    return np.nan, np.nan, np.nan, np.nan

  times = []
  for line in stdout.splitlines():
    if line.startswith("run_id=") or line.startswith("run="):
      try:
        times.append(float(line.split(",")[-1]))
      except ValueError:
        continue
    elif line.startswith("Trial"):
      try:
        times.append(float(line.split()[-1]))
      except (AttributeError, ValueError):
        continue
  if len(times) < 3:
      return np.nan, np.nan, np.nan, np.nan
  times = times[2:]  # Warmup
  return (
      geometric_mean(times) if times else np.nan,
      stdev(times) if len(times) > 1 else 0.0,
      min(times) if times else np.nan,
      max(times) if times else np.nan,
  )

def jobs_to_dataframe() -> pd.DataFrame:
  """Parses all log files in the logs/ directory and returns a DataFrame."""
  parsed_jobs = []
  log_files = glob.glob("logs/amd/*.log")

  for log_file in log_files:
    filename = os.path.basename(log_file)
    parts = filename.split(" ")
    
    implementation = parts[0]
    dataset = parts[1].removesuffix(".mtx")
    
    match = re.search(r"(\d+)cpus", parts[2])
    if not match:
      continue
    ncpus = int(match.group(1))
    chunksize = 0
    if implementation == "pthreads":
      match = re.search(r"chunk(\d+)", parts[3])
      if not match:
        continue
      chunksize = int(match.group(1))
    
    avg_runtime, std_runtime, min_runtime, max_runtime = parse_log_file(log_file)

      
    if pd.isna(avg_runtime):
      continue

    parsed_jobs.append(
      {
        "board": "amd",
        "implementation": implementation,
        "dataset": dataset,
        "num_cpus": ncpus,
        "chunksize": chunksize,
        "runtime": avg_runtime,
        "std_runtime": std_runtime,
        "min_runtime": min_runtime,
        "max_runtime": max_runtime,
      }
    )
  if not parsed_jobs:
    return pd.DataFrame()
    
  return pd.DataFrame(parsed_jobs)

if __name__ == "__main__":
  df = jobs_to_dataframe()
  if not df.empty:
    output_filename = "scripts/performance_data.csv"
    df.to_csv(output_filename, index=False)
    print(f"Data saved to {output_filename}")
  else:
    print("No log files found or processed. No data saved.")