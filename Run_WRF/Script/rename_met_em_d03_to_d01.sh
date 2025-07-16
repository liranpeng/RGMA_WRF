#!/bin/bash

# Load necessary module
conda activate nco-env
# Navigate to the target directory
cd /pscratch/sd/h/heroplr/wrf_scratch/WRF_dm/test/Domain4_20170713_chem || { echo "Directory not found"; exit 1; }

# Update metadata for each met_em.d03.* file
for f in met_em.d03.*; do
    echo "Updating metadata in $f"

    ncatted -O \
        -a grid_id,global,o,i,1 \
        -a parent_id,global,o,i,0 \
        -a i_parent_start,global,o,i,1 \
        -a j_parent_start,global,o,i,1 \
        -a i_parent_end,global,o,i,553 \
        -a j_parent_end,global,o,i,532 \
        -a parent_grid_ratio,global,o,i,1 \
        "$f"
done

# Rename all processed files from d03 to d01
echo "Renaming files from met_em.d03.* to met_em.d01.*"
for f in met_em.d03.*; do
    mv "$f" "${f/d03/d01}"
done

echo "Metadata update and renaming completed."

