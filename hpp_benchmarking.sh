#!/usr/bin/env bash
set -euo pipefail

# -------------------- Parameter-Sweeps --------------------
NP_LIST=(1 2 4 6 8 10 12)             # MPI ranks
OMP_LIST=(1 2 4 6 8 10 12)            # OpenMP threads
GRID_LIST=(64 128 256 512 1024)         # Grid sizes N
ITERS_LIST=(500 1000 2000)     # Iterations
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