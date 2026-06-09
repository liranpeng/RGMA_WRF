# WRF Low-Cloud Enhancement (MID version)

`modify_wrf_mid_lowcloud.py` perturbs WRF restart (`wrfrst`) and boundary (`wrfbdy`)
files to **increase low-cloud coverage** in the lower troposphere. It does this by
enhancing water-vapor and cloud-water mixing ratios, seeding new cloud water where
the air is nearly saturated, and re-diagnosing cloud fraction accordingly.

This is the **intermediate ("MID")** member of a perturbation ladder used for
aerosol–cloud sensitivity experiments. It sits halfway between the unmodified
control and the stronger `modify_wrf_more_lowcloud.py` perturbation:

| Parameter | Original (control) | **MID (this script)** | more_lowcloud |
|---|---|---|---|
| `QVAPOR_SCALE` | ×1.00 | **×1.10** | ×1.20 |
| `QCLOUD_SCALE` | ×1.00 | **×1.50** | ×2.00 |
| `RH_SEED_THRESHOLD` | none | **0.97** | 0.95 |
| `QCLOUD_SEED` | none | **5.0×10⁻⁶ kg/kg** | 1.0×10⁻⁵ kg/kg |

A higher seeding threshold and a smaller seed amount make MID more conservative
than `more_lowcloud`: fewer grid points get new cloud, and each gets less of it.

## What it does

The perturbation is applied only in the **low-cloud target layer** (model levels
`0`–`K_LOW-1`, roughly the surface to ~3 km). Scaling is full strength below
`K_TAPER_START` and tapers linearly back to 1.0 between `K_TAPER_START` and `K_LOW`
so there is no discontinuity at the top of the modified layer.

### Restart file (`wrfrst`)

1. **QVAPOR** — water-vapor mixing ratio scaled by `QVAPOR_SCALE` (clipped ≥ 0).
2. **QCLOUD** — cloud-water mixing ratio scaled by `QCLOUD_SCALE` in regions that
   already contain cloud (clipped ≥ 0).
3. **RH-based seeding** — relative humidity is computed from the prognostic fields
   (`T`, `P`, `PB`, `QVAPOR`). Where `RH > RH_SEED_THRESHOLD` and the existing
   `QCLOUD` is below `QCLOUD_THRESH`, new cloud water is seeded at `QCLOUD_SEED`.
4. **CLDFRA** — cloud fraction re-diagnosed within the target layer: set to 1
   wherever `QCLOUD > QCLOUD_THRESH`, and left untouched above the layer.
5. **CLDFRA_OLD** — set equal to the new `CLDFRA` inside the layer (if present).
6. **Q2** — 2-m surface moisture scaled by `QVAPOR_SCALE` (clipped ≥ 0).

### Boundary file (`wrfbdy`)

For all four lateral boundaries (`BXS`, `BXE`, `BYS`, `BYE`) and their tendency
terms (`BTXS`, `BTXE`, `BTYS`, `BTYE`):

- `QVAPOR_*` boundary values and tendencies scaled by `QVAPOR_SCALE`.
- `QCLOUD_*` boundary values and tendencies scaled by `QCLOUD_SCALE`.

Scaling the tendencies consistently keeps the boundary forcing physically coherent
with the perturbed boundary values over the integration.

## Relative-humidity calculation

RH is derived directly from WRF prognostic variables:

- Total pressure: `P_total = P + PB`
- Potential temperature: `theta = T + 300` (WRF stores `T` as `theta − 300`)
- Actual temperature: `T_act = theta · (P_total / P00)^(Rd/Cp)`
- Saturation vapor pressure: Tetens-style formula over water
- Saturation mixing ratio: `q_sat = 0.622·e_s / (P_total − 0.378·e_s)`
- `RH = QVAPOR / q_sat`

Constants used: `P00 = 1.0×10⁵ Pa`, `Rd/Cp = 0.28571428`.

## Requirements

- Python 3
- [`xarray`](https://docs.xarray.dev/)
- [`numpy`](https://numpy.org/)
- A NetCDF backend for xarray (e.g. `netCDF4`)

```bash
pip install xarray numpy netCDF4
```

## Configuration

Edit the path and parameter blocks near the top of the script before running.

**Paths**

```python
DATA_DIR   = "/pscratch/sd/h/heroplr/RGMA_homecopy/WRF_dm_Stam3/test/20170714_mid_3aer/"
WRFRST_IN  = DATA_DIR + "wrfrst_d01_2017-07-15_07:00:00"
WRFBDY_IN  = DATA_DIR + "wrfbdy_d01"
WRFRST_OUT = DATA_DIR + "wrfrst_d01_2017-07-15_07:00:00_midlowcloud"
WRFBDY_OUT = DATA_DIR + "wrfbdy_d01_midlowcloud"
```

**Layer and perturbation controls**

```python
K_LOW         = 28      # top of the low-cloud target layer (exclusive)
K_TAPER_START = 23      # full scaling below this; taper to 1.0 over [K_TAPER_START, K_LOW)

QVAPOR_SCALE  = 1.10
QCLOUD_SCALE  = 1.50

RH_SEED_THRESHOLD = 0.97     # seed only where RH exceeds this
QCLOUD_SEED       = 5.0e-6   # kg/kg of cloud water to seed
QCLOUD_THRESH     = 1.0e-7   # kg/kg; threshold for CLDFRA = 1
```

## Usage

```bash
python modify_wrf_mid_lowcloud.py
```

The script reads the input files, prints a configuration banner and a step-by-step
log (with before/after low-level means for each field), and writes the perturbed
restart and boundary files to the output paths. The inputs are **not** modified
in place — new files with the `_midlowcloud` suffix are created.

## Output

| File | Description |
|---|---|
| `wrfrst_d01_..._midlowcloud` | Perturbed restart: QVAPOR, QCLOUD, CLDFRA, CLDFRA_OLD, Q2 |
| `wrfbdy_d01_midlowcloud` | Perturbed boundaries and tendencies for QVAPOR and QCLOUD |

Point the WRF run at these files (rename or update `namelist.input` as appropriate)
to integrate the perturbed initial/boundary state.

## Notes and caveats

- The level indices (`K_LOW`, `K_TAPER_START`) are tied to the vertical grid of the
  intended case; adjust them if your domain uses a different number of model levels
  or vertical spacing.
- Cloud fraction is treated as a simple 0/1 diagnostic from `QCLOUD` within the
  target layer. This is a deliberate simplification for the perturbation experiment,
  not a physically prognostic cloud-fraction scheme.
- RH can exceed 1 in the computed field; only the saturation threshold matters for
  seeding, so this is expected and harmless here.
- The script assumes a single time record (index `0`) in the restart file.
