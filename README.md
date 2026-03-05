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
   [upstream wiki](https://github.com/TrynkaLab/tglow-pipeline/wiki/1-Installation),
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

### Slurm head-job script template

Create a script like the following and submit it with `sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=tglow
#SBATCH --account=project_XXXXXXX
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/project_XXXXXXX/my_run/tglow_%j.log

module load nextflow/25.10.2

export NXF_HOME=/scratch/project_XXXXXXX/.nextflow

# Run from the scripts directory so relative paths in your config resolve correctly
cd /scratch/project_XXXXXXX/my_run/scripts

nextflow run /projappl/project_XXXXXXX/tglow-pipeline/main.nf \
  -profile puhti \
  -c my_config.config \
  -entry run_pipeline \
  -w /scratch/project_XXXXXXX/my_run/workdir \
  -with-report logs/run_pipeline.html \
  -with-trace logs/run_pipeline.trace \
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

---

# Tglow: Nextflow pipeline for analyzing HCI data

> Check out our pre-print here: https://www.biorxiv.org/content/10.64898/2026.02.10.704860v1

This repo contains the nextflow pipeline and binaries and scripts to run a tglow-pipeline instance for the analysis of high content imaging data.
A detailed walkthrough of the steps, installation and configuration is given on the [wiki](https://github.com/TrynkaLab/tglow-pipeline/wiki) and a full list of options can be found in [docs/parameters.md](docs/parameters.md). A guided tutorial with example data is available [here](https://github.com/TrynkaLab/tglow-pipeline/wiki/7-Guided-example)

There are three components to the overall workflow
1. [tglow-pipeline](https://github.com/TrynkaLab/tglow-pipeline) - Nextflow files and Python scripts for running pipeline processes
2. [tglow-core](https://github.com/TrynkaLab/tglow-core) - A python library with IO, parsing and convenience functions based on AICSImageIO.
3. [tglow-r](https://github.com/TrynkaLab/tglow-r) - A Seurat-like R package for analyzing the output HCI features 


# Installation & dependencies
See [here](https://github.com/TrynkaLab/tglow-pipeline/wiki/1-Installation) for full install instructions of all pipeline components.

# Pipeline overview

The following readme gives a high level overview, for more detailed guide please see the wiki. The pipeline consists of two main stages:
- stage: prepare and standardize raw images into a well/field-organized OME-TIFF layout with metadata.
- run_pipeline: perform image processing and feature extraction on the staged images.

Both stages are implemented as Nextflow workflows and can be run independently using `-entry stage|run_pipeline`. 

<img src="docs/workflow.png" style="width:50%; height:auto;">


> Some steps in the pipeline require GPU's to be available. These are semgmentation and deconvolution. Deconvolution will not run without GPU. Segmentation (CellPose) will run, but we only reccomend this in cases where you are generating masks in 2D. In 3d the computational burden for large datasets will be too much for CPU. 

> The pipeline is intended to run on high performance compute (HPC) clusters, and bundled resource profiles should work for most HPC, but some tweaks to queue names and GPU settings may be required as flags differ between vendors and HPC configurations. Go to conf/processes.config and search for `queue` and `clusterOptions` to update. Furthermore each HPC is different, with different machines and resource limits. You may need to add a profile for your HPC enviroment in the conf folder. The nf-core config directory may be of help for your HPC: https://nf-co.re/configs/. If something is unclear, feel free to raise an issue on github.  

> If you dont want to run the pipeline on HPC but run it locally, supply `-profile local`.

## 1) stage
Purpose: Stage Revity/PerkinElmer (currently Phenix or Operetta) acquisitions into a reproducible plate/row/col/field.ome.tiff structure and capture metadata (channel names, pixel sizes, channel order, original index files).

-> If you don't have a Phenix or Operetta export, you can skip this step, but will need to organize the images using your own script. See more details [here](https://github.com/TrynkaLab/tglow-pipeline/wiki/3-Staging-data)

Input:
- PerkinElmer index.xml / index.idx.xml and raw instrument files (or manually organized raw files).

Output:
- plate_name/row/col/field.ome.tiff (with metadata)
- manifest listing wells to process (used as a Nextflow channel)
- auxiliary files to capture provenance (index.xml, channel maps, etc.) (optional)

Nextflow processes:
1. prepare_manifest — create a manifest with wells/fields to run (re-usable Nextflow channel)
2. fetch_raw — read raw files and write standardized OME-TIFFs and metadata

## 2) run_pipeline
Purpose: Run the core image-processing and feature-extraction steps on the staged images. The workflow is modular — many steps are optional or configurable.

Input:
- Staged OME-TIFFs (from `stage`) and optional per-plate/field metadata (flatfields, registration references, etc.)

Output:
- Segmentation outputs, registration matrices, flatfields, extracted feature tables, and logs/artifacts needed for downstream analysis.

Main processing steps (in typical execution order — each step can be enabled/disabled via config):
1. estimate flatfield (Polynomial / BaSiCPY) (optional)
   - Parallelization: per-plate + channel or single flatfield for all plates + channels.
   - Output: flatfield images only (no transformed images saved).
2. register (cross correlation / pystackreg) (optional)
   - Parallelization: per-well
   - Output: registration matrices (no transformed images saved).
3. cellpose segmentation
   - Parallelization: per-well, GPU-enabled
   - Notes: If registration is used, segmentation currently runs on the reference plate. Nucleus channel optional but segmentation is required.
   - Output: 2D or 3D cell & nucleus masks as tiffs
4. deconvolute with CLIJ2-fft (optional)
   - Parallelization: per-well, GPU-enabled
   - Output: deconvolved images (creates a data copy)
5. finalizing images
   - Parallelization: per-well
   - Applies all the registration, flatfields, scaling, max projection to the (deconvolved) images and collects the masks
   - Output: Analysis reade OME-TIFFs
6. feature extraction with CellProfiler
   - Parallelization: per-well
   - Stage images into a CellProfiler-compatible layout, apply flatfields and registration (if enabled), and run feature extraction.
   - Outputs: CellProfiler artifacts as a zip archive per well
7. cellcrops (optional)
   - Parallelization: per-well
   - Produces a HDF5 file for each field where each h5 group is a cell
   - Outputs: h5 file with fully processed cellcrops   

# Options
See nextflow.config or the [docs/parameters.md](docs/parameters.md) for available options and their descriptions.

# Quick Usage

Prerequisites:
- Nextflow and conda
- Completed [install instructions](https://github.com/TrynkaLab/tglow-pipeline/wiki/1_installation)

I strongly reccomend to configure through a configuration file, altough parameters can be overridden on the commandline. I would reccomend a project structure as follows:

- my_project
  - results: By default this is where the pipeline stores outputs
  - scripts
    - logs
    - my_config.config
    - run_pipeline.sh
  - workdir: By default this is the Nextflow workdir

Quick examples:

Stage PerkinElmer data from a raw export:
```
nextflow \
-log logs/stage.nextflow.log \
run </path/to/main.nf> \
-profile <your profile> \
-w ../workdir \
-resume \
-entry stage \
-with-report logs/stage.nextflow.html \
-with-trace logs/stage.nextflow.trace \
-c my_config.config"
```

Run the main pipeline on staged images:
```
nextflow \
-log logs/run_pipeline.nextflow.log \
run </path/to/main.nf> \
-profile <your profile> \
-w ../workdir \
-resume \
-entry run_pipeline \
-with-report logs/run_pipeline.nextflow.html \
-with-trace logs/run_pipeline.nextflow.trace \
-c my_config.config"
```

# Getting help
See [known issues and notes](https://github.com/TrynkaLab/tglow-pipeline/wiki/Kown-issues). If you find an issue please raise it on the git or contact us directly.

# Authors:
- Olivier Bakker
- Francesco Cisterno

# References
- https://github.com/clij/clij2-fft
- https://cellprofiler.org/
- https://scikit-image.org/
- https://github.com/MouseLand/cellpose
- https://github.com/glichtner/pystackreg/tree/master
- https://basicpy.readthedocs.io/en/latest/
