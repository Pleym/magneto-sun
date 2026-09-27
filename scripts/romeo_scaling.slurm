#!/bin/bash
# Phase 5 sur ROMEO : compile, teste, puis lance la campagne de mesures
# (scripts/bench_scaling.sh) sur un nœud complet.
# Avant, sur la frontale (les nœuds de calcul n'ont peut-être pas internet) :
#   scripts/fetch_gong.sh mrzqs200727t0004c2233_145.fits.gz
#   scripts/fetch_hmi.sh 2233
# Puis, depuis la racine du dépôt : sbatch scripts/romeo_scaling.slurm
#SBATCH --job-name=magnetosun-bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
# À adapter à ROMEO (noms non vérifiés ici) : partition, compte, cœurs par nœud (lscpu)
##SBATCH --partition=...
##SBATCH --account=...
set -euo pipefail

# À adapter à la sortie de « module avail » : compilateurs, CMake, BLAS/LAPACK, cfitsio.
# module load ...

cmake -S . -B build/romeo -DCMAKE_BUILD_TYPE=Release
cmake --build build/romeo -j "$SLURM_CPUS_PER_TASK"
ctest --test-dir build/romeo --output-on-failure

export OMP_PROC_BIND=close OMP_PLACES=cores
scripts/bench_scaling.sh build/romeo data/mrzqs200727t0004c2233_145.fits.gz \
  data/hmi.Synoptic_Mr.2233.fits "$SLURM_CPUS_PER_TASK" "data/bench_romeo_${SLURM_JOB_ID}.csv"
