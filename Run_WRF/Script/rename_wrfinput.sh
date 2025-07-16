#!/bin/bash
module load nco
# Define source and target directories
SRC_DIR="/scratch/07088/tg863871/WRF/test/Domain34_ndown"

# Create target directory if it doesn't exist

cd "$SRC_DIR"
cp wrfinput_d02_ndown wrfinput_d01_ndown
ncdump -v P_TOP wrfinput_d01_ndown
    ncatted -O \
        -a GRID_ID,global,o,i,1 \
        -a PARENT_ID,global,o,i,0 \
        -a I_PARENT_START,global,o,i,1 \
        -a J_PARENT_START,global,o,i,1 \
        -a PARENT_GRID_RATIO,global,o,i,1 \
        "wrfinput_d01_ndown"
cp wrfinput_d01_ndown wrfinput_d01
ncdump -v P_TOP wrfinput_d01
echo "All files copied, renamed, and metadata updated."

