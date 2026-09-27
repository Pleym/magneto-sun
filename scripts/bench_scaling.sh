#!/usr/bin/env bash
# Phase 5 : campagne de mesures de l'ajustement (bench_fit), écrite en CSV.
#   1. temps selon lmax : solveurs dense et par anneaux sur GONG, anneaux sur HMI
#   2. passage à l'échelle en threads : anneaux sur HMI (OpenMP), dense sur GONG (BLAS)
# Même script en local et sur ROMEO (via scripts/romeo_scaling.slurm).
# Usage : scripts/bench_scaling.sh <dossier de build> <carte GONG> <carte HMI> <threads max> <sortie.csv>
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage : $0 <dossier de build> <carte GONG> <carte HMI> <threads max> <sortie.csv>" >&2
  exit 2
fi
bench="$1/bench_fit"
gong="$2"
hmi="$3"
max_threads="$4"
out="$5"

# run <threads OpenMP> <threads BLAS> <carte> <solveur> <lmax> <répétitions>
run() {
  OMP_NUM_THREADS=$1 VECLIB_MAXIMUM_THREADS=$2 OPENBLAS_NUM_THREADS=$2 MKL_NUM_THREADS=$2 \
    "$bench" "$3" "$4" "$5" 90 "$6" | grep -v '^#' | sed "s|^|$(basename "$3"),$2,|" >> "$out"
}

# 1, 2, 4, ... jusqu'à max_threads (inclus)
thread_counts() {
  local t=1
  while (( t < max_threads )); do
    echo "$t"
    t=$((t * 2))
  done
  echo "$max_threads"
}

echo "# carte,threads_blas,solveur,lmax,pixels,threads_omp,etape1,etape2,etape3,etape4,total" > "$out"

# 1. Toute la machine : le dense parallélise par la BLAS, les anneaux par OpenMP
#    (BLAS séquentielle dans les blocs, sinon sur-souscription)
for l in 10 20 30 40 60 90; do run "$max_threads" "$max_threads" "$gong" dense "$l" 1; done
for l in 10 20 30 40 60 90 120 179; do run "$max_threads" 1 "$gong" rings "$l" 3; done
for l in 90 180 360 720; do run "$max_threads" 1 "$hmi" rings "$l" 1; done

# 2. Passage à l'échelle forte
for t in $(thread_counts); do run "$t" 1 "$hmi" rings 720 3; done
for t in $(thread_counts); do run "$t" "$t" "$gong" dense 60 1; done
echo "$out"
