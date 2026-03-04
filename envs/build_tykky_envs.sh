#!/bin/bash
# Run this inside an sinteractive session:
#   sinteractive --account project_2018108 --time 2:00:00 --mem 12000 --tmp 50
#
# Then:
#   bash /projappl/project_2018108/tglow-pipeline/envs/build_tykky_envs.sh

set -euo pipefail

ENVS_DIR=/projappl/project_2018108/tglow-pipeline/envs
TGLOW_ENV=/projappl/project_2018108/tglow-env
CPR_ENV=/projappl/project_2018108/cellprofiler-env

module load tykky

# Remove any partial/empty installs from previous interrupted runs
[ -d "$TGLOW_ENV" ] && [ -z "$(ls -A "$TGLOW_ENV")" ] && rmdir "$TGLOW_ENV"
[ -d "$CPR_ENV"   ] && [ -z "$(ls -A "$CPR_ENV")"   ] && rmdir "$CPR_ENV"

echo "=== Building tglow environment (Python 3.10) ==="
conda-containerize new --mamba \
    -r "$ENVS_DIR/req_tglow.txt" \
    --prefix "$TGLOW_ENV" \
    "$ENVS_DIR/env_tglow.yml"

echo "=== Building cellprofiler environment (Python 3.9) ==="
conda-containerize new --mamba \
    -r "$ENVS_DIR/req_cellprofiler.txt" \
    --prefix "$CPR_ENV" \
    "$ENVS_DIR/env_cellprofiler.yml"

echo "=== Verifying installs ==="
echo "tglow python: $("$TGLOW_ENV/bin/python" --version)"
echo "cellprofiler python: $("$CPR_ENV/bin/python" --version)"
"$CPR_ENV/bin/cellprofiler" --version 2>/dev/null | head -1 || echo "(cellprofiler version check skipped)"

echo "=== Done. Both environments built successfully. ==="
