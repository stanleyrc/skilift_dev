#!/bin/bash
# Submit with --array=0-8%2; output and worktree roots are explicit environment settings.
#SBATCH --job-name=cell_signatures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail
: "${SKILIFT_WORKTREE:?Set SKILIFT_WORKTREE to the feature checkout}"
: "${CELL_SIG_DATA:?Set CELL_SIG_DATA to the patient-data parent directory}"
: "${CELL_SIG_OUTPUT:?Set CELL_SIG_OUTPUT to a new staging parent directory}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
patients=(BWH69 BWH70 BWH75 MGH285 MGH301 MGH302 MGH303 MGH304 MGH305)
patient=${patients[$SLURM_ARRAY_TASK_ID]}
out=$CELL_SIG_OUTPUT/$patient
mkdir -p "$out"
git -C "$SKILIFT_WORKTREE" rev-parse HEAD > "$out/source_sha.txt"
/gpfs/commons/groups/imielinski_lab/Software/miniforge3/envs/mskilab_ne1/bin/Rscript "$SKILIFT_WORKTREE/inst/gos_single_cell/export_cell_signatures.R" \
  --skilift "$SKILIFT_WORKTREE" --patient_dir "$CELL_SIG_DATA/$patient" --out_dir "$out" --work_dir "$out/fit_work"
