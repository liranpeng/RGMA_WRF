# WRF / WRF-Chem Docker Environment 

This repository provides Docker and Podman configurations to build, run, and deploy **WRF** and **WRF-Chem** models on HPC systems such as NERSC Perlmutter.

The workflow is adapted from Xia Yan’s original instructions and extended to support additional configurations and output capabilities.

---

## 📂 Docker Images Overview

Four main Dockerfiles are provided:

| Dockerfile | Purpose |
|------------|---------|
| `Dockerfile_perlmutter` | Build WRF-Chem from NCAR Derecho WRF base image. |
| `Dockerfile_liran` | WRF-Chem with extra output variables. |
| `Dockerfile_nochem_mpi` | Standard WRF compiled with MPI only. |
| `Dockerfile_nochem_mpi_and_openmp` | Standard WRF compiled with MPI and OpenMP. |

---

## 🛠️ Build Instructions

**Step 1 – Prepare the Dockerfile**

In your working directory, rename or copy the desired Dockerfile to `Dockerfile`. For example:

```bash
cp Dockerfile_liran Dockerfile
Step 2 – Build the image

Use Podman to build and tag the image:

podman build -t heroplr/wrf-chem-liran:test .
This builds the image and tags it as test to avoid overwriting the latest tag.

📤 Push the Image to Docker Hub

Login to Docker Hub:

podman login docker.io
Enter your Docker Hub username and personal access token.

Push the image:

podman push heroplr/wrf-chem-liran:test
After pushing, the image is available for use on HPC clusters via Shifter.

📥 Pull the Image on HPC (Shifter)

On Perlmutter or other HPC systems supporting Shifter:

shifterimg pull docker:heroplr/wrf-chem-liran:test
You can verify that the image is available locally:

shifterimg images | grep heroplr
🖥️ Running the Container Interactively

Launch an interactive shell:

shifter --image=docker:heroplr/wrf-chem-liran:test bash
⚙️ Running WRF via SLURM Batch Job

Below is an example SLURM script to run wrf.exe using MPI and OpenMP:

#!/bin/bash
#SBATCH -N 8
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -J wrf_n8
#SBATCH -A m4334
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@domain.com
#SBATCH -L scratch,cfs
#SBATCH -C cpu
#SBATCH --image=docker:heroplr/wrf-nochem-omp

nodes_total=8
omp_threads=4
cores_per_node=128
total_mpi_ranks=$((cores_per_node * nodes_total / omp_threads))

module load cpu
module load cray-mpich

echo "Starting WRF..."

srun -N $nodes_total -n $total_mpi_ranks --mpi=cray_shasta \
     shifter --env OMP_NUM_THREADS=$omp_threads \
             --env LD_LIBRARY_PATH=/container/hdf5/lib:/container/netcdf/lib:$LD_LIBRARY_PATH \
             /container/WRF/main/wrf.exe

echo "WRF completed."
Replace YOUR_EMAIL@domain.com with your email address.

📝 Notes

Environment Variables
OMP_NUM_THREADS controls the number of OpenMP threads per MPI rank.
LD_LIBRARY_PATH is set to include container libraries.
WRF Input Files
Before running, prepare your run directory with:
namelist.input
met_em*
Any other required input files.
Image Tags
Use descriptive tags (test, latest, v1) to manage versions.
Container Paths
Binaries are located under /container/WRF/main/.
📧 Contact

For questions or contributions, please contact:

Liran Peng (liranp@uci.edu)
