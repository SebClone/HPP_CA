#!/bin/bash
#SBATCH --job-name=mpi_benchmark
#SBATCH --output=mpi_benchmark_output.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sroth@hs-koblenz.de
#SBATCH --time=0-16:00:00
#SBATCH --ntasks=10           # Maximal benötigte Prozesszahl
#SBATCH --cpus-per-task=10 # OpenMP Threads pro Prozess
#SBATCH --mem-per-cpu=2000

module purge
module load OpenMPI/5.0.7-GCC-14.2.0

#!/usr/bin/env bash
set -euo pipefail

# -------------------- Parameter-Sweeps --------------------
NP_LIST=(1 2 4 6 8 10)             # MPI ranks
OMP_LIST=(1 2 4 6 8 10)            # OpenMP threads
GRID_LIST=(32 64 128 256 1024 2048 4096)         # Grid sizes N
ITERS_LIST=(1000 5000)     # Iterations
REPEATS=1                 # number of repetitions per configuration

# -------------------- Build once ---------------------------
make clean
make BUILD=par

# -------------------- Output CSV ---------------------------
OUT_DIR="bench_results"
CSV_FILE="$OUT_DIR/hpp_bench.csv"
mkdir -p "$OUT_DIR"

# CSV-Header (einmalig)
if [[ ! -s "$CSV_FILE" ]]; then
  echo "np,omp_threads,mode,grid,iters,repeat,runtime_s" > "$CSV_FILE"
fi

# -------------------- Sweep (encrypt -> decrypt direkt) ---
for np in "${NP_LIST[@]}"; do
  for omp in "${OMP_LIST[@]}"; do
    for grid in "${GRID_LIST[@]}"; do
      for iters in "${ITERS_LIST[@]}"; do
        for rep in $(seq 1 "$REPEATS"); do
          export OMP_NUM_THREADS="$omp"

          # --- Encrypt ---
          ENC_OUT=$(make -s encrypt NP="$np" GRID="$grid" ITERS="$iters")
          enc_runtime=$(echo "$ENC_OUT" | awk '/Total runtime:/ {print $(NF-1)}' | tail -n1)
          enc_runtime=${enc_runtime:-NaN}
          echo "$np,$omp,encrypt,$grid,$iters,$rep,$enc_runtime" >> "$CSV_FILE"

          # --- Decrypt (direkt danach, gleiche Konfiguration) ---
          # Hinweis: GRID wird im Decrypt aus der Meta gelesen (GRID-Arg wird ignoriert)
          DEC_OUT=$(make -s decrypt NP="$np" ITERS="$iters")
          dec_runtime=$(echo "$DEC_OUT" | awk '/Total runtime:/ {print $(NF-1)}' | tail -n1)
          dec_runtime=${dec_runtime:-NaN}
          echo "$np,$omp,decrypt,$grid,$iters,$rep,$dec_runtime" >> "$CSV_FILE"

        done
      done
    done
  done
done

rows=$(($(wc -l < "$CSV_FILE") - 1))
echo "Wrote $rows rows to $CSV_FILE"
