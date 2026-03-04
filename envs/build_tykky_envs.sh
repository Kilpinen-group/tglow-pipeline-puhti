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

build_env() {
    local label="$1"; local prefix="$2"; shift 2
    if [ -d "$prefix" ] && [ -n "$(ls -A "$prefix")" ]; then
        echo "=== $label: already built at $prefix — skipping ==="
        return
    fi
    # Remove empty dir left by a previously interrupted build
    [ -d "$prefix" ] && rmdir "$prefix"
    echo "=== Building $label ==="
    conda-containerize new --mamba "$@" --prefix "$prefix"
}

build_env "tglow environment (Python 3.10)" "$TGLOW_ENV" \
    -r "$ENVS_DIR/req_tglow.txt" \
    "$ENVS_DIR/env_tglow.yml"

# Note: pip runs via --post-install (not -r) so the conda env is active when pip
# installs cellprofiler. This ensures pip sees conda-installed wxpython as already
# satisfied and can find mysql_config / java from conda on PATH.
build_env "cellprofiler environment (Python 3.9)" "$CPR_ENV" \
    --post-install "$ENVS_DIR/post_install_cellprofiler.sh" \
    "$ENVS_DIR/env_cellprofiler.yml"

echo "=== Verifying installs ==="
echo "tglow python: $("$TGLOW_ENV/bin/python" --version)"
echo "cellprofiler python: $("$CPR_ENV/bin/python" --version)"
"$CPR_ENV/bin/cellprofiler" --version 2>/dev/null | head -1 || echo "(cellprofiler version check skipped)"

echo "=== Done. Both environments built successfully. ==="
