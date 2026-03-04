#!/bin/bash
# Post-install script for the cellprofiler Tykky environment.
# Runs inside the container with the conda env ACTIVE, so:
#   - mysql_config is on PATH -> mysqlclient compiles
#   - java is on PATH -> python-javabridge compiles
#
# conda's wxpython is a proper binary install but does not write pip-compatible
# dist-info metadata. We create a zero-file stub wheel and install it via uv
# --find-links so the resolver sees wxPython as already satisfied and never
# tries to build it from source (which would fail — no GTK3 dev headers).
set -euo pipefail

# Install uv: modern SAT-based resolver with no backtracking depth limit.
# pip's default resolver hits resolution-too-deep on cellprofiler 4.2.8's complex tree.
pip install --quiet uv

# --- 1. Build a stub wxPython wheel (zipfile only — no setuptools needed) ---
# conda-forge's wxpython is a proper binary but lacks pip dist-info.
# We create a stub wheel with the conda version and pass it via --find-links
# so uv uses it instead of trying PyPI (where only source distributions exist).

# Get conda-installed wxpython version
WX_VER=$(python -c "
import pathlib, json, sys
fs = list(pathlib.Path(sys.prefix).glob('conda-meta/wxpython-*.json'))
if not fs:
    print('ERROR: conda wxpython not found', file=sys.stderr); sys.exit(1)
print(json.loads(fs[0].read_text())['version'])
")
echo "Conda wxPython version: ${WX_VER}"

STUB_DIR=$(mktemp -d)
STUB_WHL="${STUB_DIR}/wxPython-${WX_VER}-py3-none-any.whl"

python << PYEOF
import zipfile, pathlib
ver, whl_path = "${WX_VER}", "${STUB_WHL}"
di = f"wxPython-{ver}.dist-info"
with zipfile.ZipFile(whl_path, "w") as z:
    z.writestr(f"{di}/METADATA",
        f"Metadata-Version: 2.1\nName: wxPython\nVersion: {ver}\n")
    z.writestr(f"{di}/WHEEL",
        "Wheel-Version: 1.0\nGenerator: stub\nRoot-Is-Purelib: true\nTag: py3-none-any\n")
    z.writestr(f"{di}/INSTALLER", "pip\n")
    z.writestr(f"{di}/RECORD",
        f"{di}/METADATA,,\n{di}/WHEEL,,\n{di}/INSTALLER,,\n{di}/RECORD,,\n")
print(f"Stub wheel written: {whl_path}")
PYEOF

# --- 2. Build constraint file ---
# Pin wxPython to the conda-installed version so uv uses the local stub wheel
# rather than attempting a PyPI source build.
# Everything else is resolved freely by uv (see overrides in step 3).
CONSTRAINTS="${STUB_DIR}/constraints.txt"
echo "wxPython==${WX_VER}" > "${CONSTRAINTS}"
echo "  wxPython==${WX_VER} (constraint)"

# --- 3. Install cellprofiler + tglow-core using uv ---
# cellprofiler==4.2.8 pins deps that conflict with tglow-core's requirements:
#   scikit-image: cellprofiler==0.18.3  vs  tglow-core>=0.20.0
#   numpy:        cellprofiler<1.25     vs  tglow-core>=1.26.4
# The upstream wiki confirms cellprofiler works fine with newer versions at
# runtime. We use uv --override to make the resolver accept the newer versions.
# numpy is capped at <2.0 because python-javabridge uses the pre-2.0 numpy C API.
OVERRIDES="${STUB_DIR}/overrides.txt"
printf 'scikit-image>=0.20.0\nnumpy>=1.26.4,<2.0\n' > "${OVERRIDES}"

# --system:             install into the active conda env (no virtualenv required)
# --find-links:         local stub wheel for wxPython
# --constraint:         pin wxPython to the conda/stub version
# --override:           accept newer numpy/scikit-image than cellprofiler pins
# --no-build-isolation: use conda env's numpy/cython when building C extensions
#                       (mysqlclient needs mysql_config; python-javabridge needs java)
uv pip install --system --no-build-isolation \
    --find-links "${STUB_DIR}" \
    --constraint "${CONSTRAINTS}" \
    --override "${OVERRIDES}" \
    tglow-core cellprofiler==4.2.8

rm -rf "${STUB_DIR}"
