# CSC Puhti setup for tglow-pipeline

> This section documents Puhti-specific installation of the [tglow-pipeline](https://github.com/TrynkaLab/tglow-pipeline) developed by the Trynka lab. We refer to the TrynkaLab tglow-pipeline [wiki](https://github.com/TrynkaLab/tglow-pipeline/wiki) for full instructions on how to install, stage, configure and run the tglow-pipeline.

## Prerequisites

- Access to CSC Puhti with a project account (e.g. `project_XXXXXXX`)
- GPU partition access (required for cellpose segmentation and deconvolution)
- Nextflow: `module load nextflow/25.10.2` (or the latest available version)

## Directory layout

We recommend the following layout, substituting your project number:

```
/projappl/project_XXXXXXX/
├── tglow-pipeline/          # this repository
├── tglow-env/               # Tykky env: Python 3.10 (tglow-core, cellpose, …)
└── cellprofiler-env/        # Tykky env: Python 3.9 (cellprofiler 4.2.8)

/scratch/project_XXXXXXX/
└── my_run/
    ├── pipeline_testdata/   # input data
    ├── workdir/             # Nextflow work directory
    └── results/             # pipeline outputs
```

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-org>/tglow-pipeline.git \
    /projappl/project_XXXXXXX/tglow-pipeline
```

### 2. Build the Tykky environments

Puhti requires software to be containerized via
[Tykky](https://docs.csc.fi/computing/containers/tykky/) rather than installed
as plain conda environments (Lustre filesystem limitation). A build script is
provided:

```bash
# Edit the paths at the top of the script to match your project number
nano /projappl/project_XXXXXXX/tglow-pipeline/envs/build_tykky_envs.sh

# Start an interactive compute session and run the script
sinteractive --account project_XXXXXXX --time 1:30:00 --mem 16000 --tmp 100
bash /projappl/project_XXXXXXX/tglow-pipeline/envs/build_tykky_envs.sh
```

This builds two environments sequentially (~60–90 min total):

| Environment | Python | Key packages |
|-------------|--------|--------------|
| `tglow-env` | 3.10 | tglow-core, cellpose 3.0.8, clij2-fft, RedLionfish, h5py |
| `cellprofiler-env` | 3.9 | cellprofiler 4.2.8, tglow-core |

**Important — cellprofiler-env build notes:**

`cellprofiler-env` uses `--post-install` (not `-r`) so the post-install script
runs with the full conda environment active. The script (`envs/post_install_cellprofiler.sh`)
does the following:

1. **wxPython stub wheel.** `wxPython` has no pre-built Linux wheels on PyPI
   and cannot be built from source inside the container (no GTK3 dev headers).
   conda-forge provides a working binary. The post-install script creates a
   zero-file stub `.whl` using Python's `zipfile` module (no setuptools needed)
   and installs it with `--find-links` so the package manager sees wxPython as
   already satisfied at the conda version. A constraint file pins it to that
   exact version to prevent upgrades.

2. **System tools for C-extension builds.** `mysql_config` (from conda `mysql`)
   and `java` (from conda `openjdk`) must be on PATH when compiling `mysqlclient`
   and `python-javabridge`. The `--no-build-isolation` flag ensures these conda
   packages are visible during the build.

3. **uv resolver.** `cellprofiler==4.2.8` has a dependency tree too complex for
   pip's backtracking resolver (hits the depth limit). The post-install script
   installs `uv` first and uses `uv pip install --system` instead, which uses a
   SAT-based resolver without a depth limit.

4. **Dependency overrides.** `cellprofiler==4.2.8` pins `scikit-image==0.18.3`
   and requires `numpy<1.25`, but `tglow-core` needs `scikit-image>=0.20.0` and
   `numpy>=1.26.4`. Following the approach in the
   [tglow wiki](https://github.com/TrynkaLab/tglow-pipeline/wiki/1-Installation),
   cellprofiler works fine with newer versions at runtime; `uv --override` tells
   the resolver to accept `scikit-image>=0.20.0` and `numpy>=1.26.4,<2.0`.
   (numpy must stay below 2.0 because `python-javabridge 4.0.4` uses the old
   numpy C API.)

**Note:** The `tglow-core` package installs as the Python module `tglow`
(i.e. `import tglow`, not `import tglow_core`).

After the build completes, verify both environments:

```bash
# cellprofiler-env: check cellprofiler, tglow, and wx are importable
/projappl/project_XXXXXXX/cellprofiler-env/bin/python -c "
import cellprofiler, tglow, wx
print('cellprofiler:', cellprofiler.__version__)
print('tglow:        ok')
print('wx:          ', wx.__version__)
"

# tglow-env
/projappl/project_XXXXXXX/tglow-env/bin/python -c "import tglow; print('tglow: ok')"
```

A scipy UserWarning about numpy version on import is expected and harmless.

### 3. Configure the Puhti profile

The Puhti Nextflow profile is in `conf/puhti.config`. Before using it, update
the two path variables at the top of the file to match your project:

```groovy
// conf/puhti.config
TGLOW_ENV = '/projappl/project_XXXXXXX/tglow-env'
CPR_ENV   = '/projappl/project_XXXXXXX/cellprofiler-env'
```

Also update the billing account in the `clusterOptions` lines:

```groovy
clusterOptions = '--account=project_XXXXXXX'
// and for GPU processes:
clusterOptions = '--account=project_XXXXXXX --gres=gpu:v100:1'
```

The profile maps tglow's internal resource labels to Puhti partitions:
- CPU labels (`small`, `normal`, `himem`, `*_img`) → `small` partition
- GPU labels (`gpu_*`) → `gpu` partition with one V100

### 4. Register the profile in nextflow.config

If it is not already present, add the Puhti profile inside the `profiles {}`
block in `nextflow.config`:

```groovy
profiles {
    // ... existing profiles ...
    puhti { includeConfig 'conf/puhti.config' }
}
```

## Running a pipeline

### Slurm head-job script

We recommend organising each experiment as follows, then submitting with
`sbatch run_<experiment>.sh`:

```
/scratch/project_XXXXXXX/my_experiment/
├── workdir/              # Nextflow work directory (-w flag)
├── results/              # pipeline outputs (rn_publish_dir / rn_image_dir)
├── data/                 # raw instrument data (stage entry only)
└── scripts/
    ├── logs/             # Nextflow logs, HTML reports, trace files
    └── inputs/
        ├── manifest.tsv          # plate/well/field manifest
        ├── control_list.tsv      # control wells for rn_autoscale (if used)
        └── my_experiment.config  # per-run Nextflow params
```

The run script below is self-documenting — edit the `UPDATE FOR EACH RUN`
sections and leave the rest unchanged:

```bash
#!/bin/bash

# ===========================================================================
# tglow-pipeline run script for CSC Puhti
#
# To start a new run:
#   1. Copy this script and the scripts/ directory to your experiment folder
#   2. Edit the "UPDATE FOR EACH RUN" sections below
#   3. Edit inputs/my_experiment.config (manifests, output dirs, pipeline params)
#   4. Submit: sbatch run_my_experiment.sh
# ===========================================================================

# --- UPDATE FOR EACH RUN: Slurm job identity ---
#SBATCH --job-name=my_experiment
#SBATCH --account=project_XXXXXXX
#SBATCH --output=/scratch/project_XXXXXXX/my_experiment/scripts/logs/my_experiment_%j.log
# ------------------------------------------------

# Fixed Slurm settings (head job only — child jobs are submitted by Nextflow)
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,FAIL,END

# ===========================================================================
# UPDATE FOR EACH RUN
# ===========================================================================

# Pipeline entry point: "stage" (raw → OME-TIFF) or "run_pipeline" (process staged images)
ENTRY="run_pipeline"

# Short name for this run — used in report/trace filenames
RUN_NAME="my_experiment_run"

# Directory containing inputs/ and logs/ for this run.
# Nextflow is launched from here so relative paths in CONFIG resolve correctly.
SCRIPTS_DIR="/scratch/project_XXXXXXX/my_experiment/scripts"

# Per-run Nextflow config: manifests, output dirs, and pipeline parameters.
CONFIG="${SCRIPTS_DIR}/inputs/my_experiment.config"

# Nextflow work directory for intermediate files.
# Shared across runs when using -resume; use a separate dir to start fresh.
WORKDIR="/scratch/project_XXXXXXX/my_experiment/workdir"

# ===========================================================================

# Pipeline settings (update if pipeline installed in a different location)
PIPELINE="/projappl/project_XXXXXXX/tglow-pipeline/main.nf"

module load nextflow/25.10.2

export NXF_HOME=/scratch/project_XXXXXXX/.nextflow

mkdir -p "${SCRIPTS_DIR}/logs"
cd "${SCRIPTS_DIR}"

nextflow run "${PIPELINE}" \
  -profile puhti \
  -c "${CONFIG}" \
  -entry "${ENTRY}" \
  -w "${WORKDIR}" \
  -with-report "logs/${RUN_NAME}.nextflow.html" \
  -with-trace "logs/${RUN_NAME}.nextflow.trace" \
  -resume
```

Use `-entry stage` instead of `-entry run_pipeline` if you first need to stage
raw PerkinElmer instrument data.

### Demo run

A ready-to-use demo script is provided at
`/scratch/project_XXXXXXX/tglow_example/run_tglow_demo.sh`. Download the demo
data first:

```bash
mkdir -p /scratch/project_XXXXXXX/tglow_example
cd /scratch/project_XXXXXXX/tglow_example
wget https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/652/S-BSST2652/Files/TEST_DATA/pipeline_testdata_v1.zip
unzip pipeline_testdata_v1.zip
sbatch run_tglow_demo.sh
```

## Puhti-specific notes and known issues

### Local scratch (`rn_scratch`)

**Do not use `rn_scratch = true` without also requesting NVMe storage from
SLURM.** The `finalize` and `cellprofiler` processes write large temporary
files (several GB per well), and without an explicit `--gres=nvme:N` allocation
the compute node only has a small default `/tmp`, which fills up immediately.

The safest option for most runs is:

```groovy
// in your experiment .config
rn_scratch = false   // write directly to Lustre scratch
```

If you want to use local NVMe for performance, add the required SLURM flag to
the relevant process labels in `conf/puhti.config`:

```groovy
withLabel: normal { clusterOptions = '--account=project_XXXXXXX --gres=nvme:50' }
```

and set `rn_scratch = true` in your run config.

### scikit-image API compatibility (`selem` → `footprint`)

CellProfiler 4.2.8 uses the old scikit-image `selem` keyword argument in its
`MeasureGranularity` module, which was renamed to `footprint` in
scikit-image ≥ 0.20. Since `tglow-core` requires scikit-image ≥ 0.20, the two
packages conflict at runtime.

**Workaround:** replace the `cellprofiler` entry-point symlink in
`cellprofiler-env/_bin/` with a small Python wrapper that monkey-patches
`skimage.morphology` before launching CellProfiler:

```bash
# Run once after building the cellprofiler-env
rm /projappl/project_XXXXXXX/cellprofiler-env/_bin/cellprofiler
```

Then create `/projappl/project_XXXXXXX/cellprofiler-env/_bin/cellprofiler`
with the following content (note: the shebang is valid only inside the
Singularity container that Tykky creates; save with **LF line endings**):

```python
#!/PUHTI_TYKKY_SYR94yI/miniforge/envs/env1/bin/python3
import sys, skimage.morphology

_orig_erosion = skimage.morphology.erosion
def _erosion(image, footprint=None, selem=None, **kwargs):
    if selem is not None and footprint is None:
        footprint = selem
    return _orig_erosion(image, footprint=footprint, **kwargs)
skimage.morphology.erosion = _erosion

_orig_dilation = skimage.morphology.dilation
def _dilation(image, footprint=None, selem=None, **kwargs):
    if selem is not None and footprint is None:
        footprint = selem
    return _orig_dilation(image, footprint=footprint, **kwargs)
skimage.morphology.dilation = _dilation

_orig_reconstruction = skimage.morphology.reconstruction
def _reconstruction(seed, mask, method='dilation', footprint=None, selem=None, offset=None):
    if selem is not None and footprint is None:
        footprint = selem
    return _orig_reconstruction(seed, mask, method=method, footprint=footprint, offset=offset)
skimage.morphology.reconstruction = _reconstruction

from cellprofiler.__main__ import main
sys.exit(main())
```

```bash
chmod +x /projappl/project_XXXXXXX/cellprofiler-env/_bin/cellprofiler
# Verify no CRLF line endings (critical — CRLF causes "bad interpreter" error):
file /projappl/project_XXXXXXX/cellprofiler-env/_bin/cellprofiler
# Should print: "Python script, ASCII text executable"
```

**Why `_bin/` and not `bin/`?** Tykky bind-mounts `_bin/` over `bin/` inside
the Singularity container (`common.sh` line 54), so scripts placed in `_bin/`
are what actually runs. The `_bin/` symlinks point to
`/PUHTI_TYKKY_SYR94yI/miniforge/...` which resolves to the squashfs only
inside the container.

---

For full documentation on pipeline configuration, parameters, and usage see the
[tglow-pipeline wiki](https://github.com/TrynkaLab/tglow-pipeline/wiki) and the
[pre-print](https://www.biorxiv.org/content/10.64898/2026.02.10.704860v1).
