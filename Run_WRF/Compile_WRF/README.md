WRF / WRF-Chem Docker Environment

This repository provides Docker and Podman configurations to build, run, and deploy WRF and WRF-Chem models on HPC systems such as NERSC Perlmutter.

The workflow is adapted from Xia Yan’s original instructions and extended to support additional configurations and output capabilities.

----------------------------------------
Docker Images Overview

Four main Dockerfiles are provided:

Dockerfile                          | Purpose
------------------------------------|---------------------------------------------------
Dockerfile_perlmutter               | Build WRF-Chem from NCAR Derecho WRF base image.
Dockerfile_liran                    | WRF-Chem with extra output variables.
Dockerfile_nochem_mpi               | Standard WRF compiled with MPI only.
Dockerfile_nochem_mpi_and_openmp    | Standard WRF compiled with MPI and OpenMP.

----------------------------------------
Build Instructions

Step 1 – Prepare the Dockerfile

In your working directory, rename or copy the desired Dockerfile to Dockerfile. For example:

cp Dockerfile_nochem_mpi_and_openmp Dockerfile

Step 2 – Build the Image

Use Podman to build and tag the image:

podman build -t heroplr/wrf-nochem-omp:test .

This builds the image and tags it as "test" to avoid overwriting "latest".

----------------------------------------
Push the Image to Docker Hub

Login to Docker Hub:

podman login docker.io

Enter your Docker Hub username and personal access token.

Push the image:

podman push heroplr/wrf-nochem-omp:test

After pushing, the image is available for use on HPC clusters via Shifter.

----------------------------------------
Pull the Image on HPC (Shifter)

On Perlmutter or other HPC systems supporting Shifter:

shifterimg pull docker:heroplr/wrf-nochem-omp:test

Verify the image:

shifterimg images | grep heroplr

----------------------------------------
Running the Container Interactively

Launch an interactive shell:

shifter --image=docker:heroplr/wrf-nochem-omp:test bash

----------------------------------------
Running WRF via SLURM Batch Job

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
#SBATCH --image=docker:heroplr/wrf-nochem-omp:test

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

----------------------------------------
Benchmark Performance

The following table summarizes performance testing conducted on NERSC Perlmutter HPC.
Note: These speed tests were based on the image built with Dockerfile_nochem_mpi_and_openmp.

Experiment | Nodes x OMP Threads | Model Minutes Simulated | Wall Clock Minutes | Model min / Real min | Hours to Simulate 1 Day | Node-Hours to Simulate 1 Day
---------- |----------------------|-------------------------|---------------------|-----------------------|--------------------------|------------------------------
1          | 4 x 1                | 70                      | 30                  | 2.33                  | 10.3                    | 41.2
2          | 4 x 2                | 90                      | 30                  | 3.0                   | 8.0                     | 32.0
3          | 4 x 4                | 100                     | 30                  | 3.33                  | 7.2                     | 28.8
4          | 4 x 8                | 110                     | 30                  | 3.67                  | 6.55                    | 26.2
5          | 8 x 1                | 80                      | 30                  | 2.67                  | 9.0                     | 72.0
6          | 8 x 2                | 120                     | 30                  | 4.0                   | 6.0                     | 48.0
7          | 8 x 4                | 120                     | 30                  | 4.0                   | 6.0                     | 48.0
8          | 8 x 8                | 120                     | 30                  | 4.0                   | 6.0                     | 48.0
9          | 12 x 2               | 105                     | 30                  | 3.5                   | 6.86                    | 82.3
10         | 12 x 4               | 120                     | 30                  | 4.0                   | 6.0                     | 72.0
11         | 12 x 8               | 128                     | 30                  | 4.27                  | 5.63                    | 67.6
12         | 16 x 2               | 120                     | 30                  | 4.0                   | 6.0                     | 96.0
13         | 16 x 4               | 120                     | 30                  | 4.0                   | 6.0                     | 96.0
14         | 16 x 8               | 132                     | 30                  | 4.4                   | 5.45                    | 87.2

----------------------------------------
Notes

Environment Variables
- OMP_NUM_THREADS controls the number of OpenMP threads per MPI rank.
- LD_LIBRARY_PATH is set to include container libraries.

WRF Input Files
Before running, prepare your run directory with:
- namelist.input
- met_em*
- Any other required input files.

Image Tags
Use descriptive tags (test, latest, v1) to manage versions.

Container Paths
Binaries are located under:
/container/WRF/main/

----------------------------------------
Contact

For questions or contributions, please contact:

Liran Peng
liranp@uci.edu
