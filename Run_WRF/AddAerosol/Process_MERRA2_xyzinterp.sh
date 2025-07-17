#!/bin/bash

# === User Settings ===
year=2017
month=6
mon=$(printf "%02d" $month)

# === Paths ===
obs_path="/pscratch/sd/h/heroplr/wrf_scratch/MERRA2_Data/"
wrf_path="/pscratch/sd/h/heroplr/wrf_scratch/MERRA2_Data/Liran_test/"
output_dir="/pscratch/sd/h/heroplr/wrf_scratch/MERRA2_Data/domain4"
mkdir -p "$output_dir"

# === Slurm job script path ===
job_script="$output_dir/submit_merra2_monthly.sh"

# === Write Slurm header ===
cat << EOF > "$job_script"
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH -J merra2_v2
#SBATCH -A m4334
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=liranp@uci.edu
#SBATCH -L scratch,cfs
#SBATCH -C cpu

conda activate ncl_stable

EOF

# === Get number of days in the month ===
days_in_month=$(cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}')

# === Loop through each day and time step ===
for day in $(seq -w 1 $days_in_month); do
    for t in {0..7}; do

        obs_file="MERRA2_400.inst3_3d_aer_Nv.${year}${mon}${day}.nc4"
        ncl_script="$output_dir/process_merra2_wrf_log_${year}${mon}${day}_t${t}.ncl"
        log_file="$output_dir/log_${year}${mon}${day}_log_t${t}.txt"

        # === Generate the NCL script ===
        cat << EOF_NCL > "$ncl_script"
load "\$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "\$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
    obs_path   = "$obs_path"
    obs_file   = "$obs_file"
    wrf_path   = "$wrf_path"
    wrfi_file  = "wrfinput_d04"
    wrf_file   = "met_em.d04.2016-07-01_08:50:00.nc"
    wrfb_file  = "wrfbdy_d04"

    f_obs = addfile(obs_path + obs_file, "r")
    olon  = f_obs->lon
    olat  = f_obs->lat
    dp    = f_obs->DELP($t, :, :, :)
    Ps    = f_obs->PS($t, :, :)

    oBC   = f_obs->BCPHILIC($t, :, :, :) + f_obs->BCPHOBIC($t, :, :, :)
    oOC   = f_obs->OCPHILIC($t, :, :, :) + f_obs->OCPHOBIC($t, :, :, :)
    oSO2  = f_obs->SO2($t, :, :, :)
    oSO4  = f_obs->SO4($t, :, :, :)

    f_wrfi = addfile(wrf_path + wrfi_file, "r")
    P_wrf = f_wrfi->P(0, :, :, :) + f_wrfi->PB(0, :, :, :)
    nlat_wrf = dimsizes(P_wrf(0,:,0))
    nlon_wrf = dimsizes(P_wrf(0,0,:))
    lev_wrf  = dimsizes(P_wrf(:,0,0))

    f_wrf = addfile(wrf_path + wrf_file, "r")
    LON_wrf = f_wrf->XLONG_M(0, :, :)
    LAT_wrf = f_wrf->XLAT_M(0, :, :)

    ;--- Get the number of vertical levels
    levs = dimsizes(dp)
    nlev = levs(0)
    nlat = levs(1)
    nlon = levs(2)

    ;--- Loop through levels to compute cumulative sum from top to bottom
    dp_cumsum = dp
    dp_cumsum(0,:,:) = dp(0,:,:) / 2.0
    do k = 1, nlev-1
      dp_cumsum(k,:,:) = dp_cumsum(k-1,:,:) + 0.5 * (dp(k,:,:)+dp(k-1,:,:))
    end do

    ; Allocate output
    pressure_mid = dp    ; same shape as delp
    pressure_mid = 0.0

    ; Compute pressure edges and mid-layers
    ptop = 2.0             ; Pa

    ; Loop over time
     p_edge = new((/nlev+1,nlat,nlon/), "float")
     p_edge(0,:,:) = ptop
     do k=1,nlev
        p_edge(k,:,:) = p_edge(k-1,:,:) + dp(k-1,:,:)
     end do
     PLm = new((/nlev, nlat, nlon/), float)
     do k=0,nlev-1
        PLm(k,:,:) = 0.5*(p_edge(k,:,:) + p_edge(k+1,:,:))
     end do

    nBC_interp  = new((/nlev, nlat_wrf, nlon_wrf/), float)
    nOC_interp  = new((/nlev, nlat_wrf, nlon_wrf/), float)
    nSO2_interp = new((/nlev, nlat_wrf, nlon_wrf/), float)
    nSO4_interp = new((/nlev, nlat_wrf, nlon_wrf/), float)
    PLm_interp  = new((/nlev, nlat_wrf, nlon_wrf/), float)
    
    do k = 0, nlev-1
        kr = nlev-1-k
        nBC_interp(kr, :, :)  = rgrid2rcm(olat, olon, oBC(k, :, :), LAT_wrf, LON_wrf, 0)
        nOC_interp(kr, :, :)  = rgrid2rcm(olat, olon, oOC(k, :, :), LAT_wrf, LON_wrf, 0)
        nSO2_interp(kr, :, :) = rgrid2rcm(olat, olon, oSO2(k, :, :), LAT_wrf, LON_wrf, 0)
        nSO4_interp(kr, :, :) = rgrid2rcm(olat, olon, oSO4(k, :, :), LAT_wrf, LON_wrf, 0)
        PLm_interp(kr, :, :)  = rgrid2rcm(olat, olon, PLm(k, :, :), LAT_wrf, LON_wrf, 0)
    end do

    linlog = 2

    ; Reverse order to go from top down
    nBC = int2p_n(PLm_interp, nBC_interp, P_wrf, linlog, 0)

    ; You can do the same for the other variables
    nOC = int2p_n(PLm_interp, nOC_interp, P_wrf, linlog, 0)

    nSO2 = int2p_n(PLm_interp, nSO2_interp, P_wrf, linlog, 0)

    nSO4 = int2p_n(PLm_interp, nSO4_interp, P_wrf, linlog, 0)

    out_file = "${output_dir}/regridded_merra2_d04_to_wrf_logv3_${year}${mon}${day}_t${t}.nc"
    f_out = addfile(out_file, "c")
    filedimdef(f_out, (/"lev", "lat", "lon"/), (/lev_wrf, nlat_wrf, nlon_wrf/), (/False, False, False/))

    f_out->lon  = LON_wrf
    f_out->lat  = LAT_wrf
    f_out->lev  = P_wrf
    f_out->nBC_interp  = nBC
    f_out->nOC_interp  = nOC
    f_out->nSO2_interp = nSO2
    f_out->nSO4_interp = nSO4

    print("Time Step $t for Day ${year}${mon}${day} Complete!")

end
EOF_NCL

        # === Append NCL call to job script ===
        echo "ncl \"$ncl_script\" > \"$log_file\" 2>&1 &" >> "$job_script"

    done
done

# === Add wait and final message ===
echo "wait" >> "$job_script"
echo "echo \"All NCL scripts completed.\"" >> "$job_script"

# === Make executable and submit ===
chmod +x "$job_script"
echo "Submitting job script: $job_script"
sbatch "$job_script"



