#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"

controller_job_name="${CONTROLLER_JOB_NAME:-shortread_snake}"
controller_partition="${CONTROLLER_PARTITION:-c_medium}"
controller_cpus="${CONTROLLER_CPUS_PER_TASK:-1}"
controller_mem="${CONTROLLER_MEM:-8G}"
controller_time="${CONTROLLER_TIME:-2-00:00:00}"
controller_log_dir="${CONTROLLER_LOG_DIR:-${repo_root}/shortread/snake/logs/controller}"
conda_env_name="${CONDA_ENV_NAME:-spliceosome-shortread}"
snakemake_jobs="${SNAKEMAKE_JOBS:-50}"
snakemake_cores="${SNAKEMAKE_CORES:-400}"
snakemake_latency_wait="${SNAKEMAKE_LATENCY_WAIT:-60}"

mkdir -p "${controller_log_dir}"

snakemake_args=(
  -s shortread/snake/Snakefile
  --jobs "${snakemake_jobs}"
  --cores "${snakemake_cores}"
  --latency-wait "${snakemake_latency_wait}"
  --executor cluster-generic
  --cluster-generic-submit-cmd "sbatch {params.cluster}"
  --rerun-incomplete
  --rerun-triggers mtime input params
  "$@"
)

printf -v repo_root_q "%q" "${repo_root}"
printf -v snakemake_cmd "%q " snakemake "${snakemake_args[@]}"

wrap_script=$(cat <<EOF
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${conda_env_name}
cd ${repo_root_q}
${snakemake_cmd}
EOF
)

printf -v wrap_q "%q" "${wrap_script}"

sbatch \
  --job-name "${controller_job_name}" \
  --qos "${controller_partition}" \
  --cpus-per-task "${controller_cpus}" \
  --mem "${controller_mem}" \
  --time "${controller_time}" \
  --output "${controller_log_dir}/${controller_job_name}.%j.log" \
  --chdir "${repo_root}" \
  --wrap "bash -lc ${wrap_q}"
