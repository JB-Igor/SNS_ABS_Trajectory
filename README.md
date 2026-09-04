# SNS ABS Trajectory

Last reviewed: 2026-09-04

`SNS_ABS_Trajectory` is an Igor Pro 9 package for trajectory-resolved
Andreev-bound-state (ABS) and local-density-of-states (LDOS) calculations in
ballistic superconductor-normal-superconductor (SNS) junctions. Real-space
masks define the normal-region geometry, from which the package constructs
semiclassical channels, solves their ABS spectra, and assembles point, line,
field-dependent, and spatially resolved LDOS results.

The repository includes a synthetic, data-free demonstration that runs the
complete workflow from geometry construction through calculated waves and
display windows.

## Features

- Channel construction from scaled two-dimensional masks.
- Two-dimensional and finite-height three-dimensional ray geometries.
- Ballistic SNS ABS spectra with finite interface transparency.
- Orbital magnetic phase and optional vortex phase contributions.
- Raw and experimentally broadened LDOS calculations.
- `LDOS(E,B)` at selected positions and along spatial lines.
- Full raw and convolved `LDOS(x,y,E)` cubes at fixed magnetic field.
- Ray, path-length-distribution, and LDOS display helpers.
- Native Igor display fallback when SIDAM is not installed.

## Physical scope

Each channel is represented by a normal-region path length, magnetic lever
arm, transparency, statistical weight, and boundary intersections. The solver
maps these channels onto quasi-one-dimensional ballistic SNS quantization
conditions and accumulates their spectral contributions.

This is a semiclassical trajectory model, not a general microscopic BdG or
diffusive Usadel solver. In particular, the three-dimensional option extends
the ray geometry through a finite junction height; it does not solve a full
three-dimensional BdG Hamiltonian.

## Requirements

- [Igor Pro 9](https://www.wavemetrics.com/software/igor-pro-900).
- No external data files are required for the included demo.

[SIDAM](https://github.com/yuksk/SIDAM) is optional. If available, the display
helpers use its color tables and presentation utilities. If it is absent, the
package compiles and the demo runs using native Igor display commands and
color tables. SIDAM 9.8.8 is the version used with the package during
development.

## Installation

1. Download or clone this repository.
2. Place the `SNS_ABS_Trajectory` directory under your Igor Pro 9 User Files
   `User Procedures` directory. In Igor, **Help > Show Igor Pro User Files**
   opens the corresponding user-files location.
3. Add the loader to the current experiment's procedure window:

   ```igorpro
   #include "SNS_ABS_Trajectory Loader"
   ```

4. Compile the procedures.
5. Choose **Macros > Load SNS_ABS_Trajectory**.

The same menu provides **Unload SNS_ABS_Trajectory**. The optional legacy
S-S'-S procedure is not loaded by default.

## Run the included demo

After loading the package, open
`Demo/SNS_ABS_Demo_ModelGeometryAndLDOS.txt` as an Igor notebook and execute
its sections from top to bottom.

On the first run in an experiment, select the distributed `config` directory
when prompted. The notebook creates an Igor named path and reuses it on later
runs in that experiment. The settings are read from
`config/SNS_DefaultSettings.tsv`.

The demo constructs a rectangular SNS junction and calculates:

1. The model mask and selected center/edge positions.
2. Central-position rays and their normal-path-length distribution.
3. Two-dimensional `LDOS(E,B)` at the center and near an edge.
4. Finite-height three-dimensional rays and `LDOS(E,B)` at both positions.
5. Top and side projections of representative three-dimensional rays.
6. Zero-field two-dimensional line DOS across the junction.
7. Raw and convolved two-dimensional surface-state `LDOS(x,y,E)` cubes for a
   centered free vortex on a 9 by 9 spatial grid.

Results are stored below:

```text
root:SNS_ABS_Demo_ModelGeometryAndLDOS
```

The notebook retains raw and convolved scientific waves separately and creates
normalized copies only for display. Re-running it intentionally deletes and
rebuilds this demo data folder and its named display windows.

## Settings and units

Public geometry-facing interfaces use:

- lengths and coordinates in nanometres (`*_nm`);
- energies in electronvolts (`*_eV`);
- magnetic fields in tesla (`*_T`);
- angles in radians (`*_rad`) or degrees (`*_deg`) as indicated by names.

The distributed TSV contains documented starting values, including the energy
and field grids, superconducting gap, Fermi energy, interface barrier,
temperature, modulation broadening, and vortex controls. Review these values
before applying the model to a different material or geometry.

Canonical channel waves include:

```text
L_N_List_nm
W_eff_List_nm
T_eff_List
wChan
Hit1x_List_nm, Hit1y_List_nm
Hit2x_List_nm, Hit2y_List_nm
```

`wChan` is a relative geometrical/statistical channel weight normalized by the
DOS construction; it is not an independently calibrated absolute DOS scale.

## Principal workflow functions

The demo is the maintained executable example. Its main public calls include:

```igorpro
SNS_LoadSettingsConfig(...)
SNS_InitDefaultSettings(...)
SNS_ExtractModesForFolder(...)
SNS_ExtractModes3DForFolder(...)
SNS_ComputeDOS_FromSettings(...)
SNS_ApplyDOS_Broadening_TplusMod(...)
SNS_LDOSmap_BFixed_FromMask(...)
SNS_MakeLDOS_E_MapFromSpectra(...)
SNS_DisplayWithScales(...)
```

For a fixed-field spatial calculation,
`SNS_LDOSmap_BFixed_FromMask(...)` produces the point coordinates, validity
mask, integrated LDOS, and complete raw/convolved local spectra.
`SNS_MakeLDOS_E_MapFromSpectra(...)` then maps those point spectra onto an
`x` by `y` by energy image stack. Preserve both the raw and convolved cubes
when comparing model physics with experimental resolution.

## Repository layout

```text
SNS_ABS_Trajectory/
  config/
    SNS_DefaultSettings.tsv
  Demo/
    SNS_ABS_Demo_ModelGeometryAndLDOS.txt
  procedures/
    SNS_ABS_Trajectory Loader.ipf
    SNS_Core.ipf
    SNS_GeometryFromMask.ipf
    SNS_Solver.ipf
    SNS_DOS.ipf
    SNS_SpatialMaps.ipf
    SNS_Broadening.ipf
    SNS_DisplayHelpers.ipf
    ...
  LICENSE
  README.md
```

## License

Released under the [MIT License](LICENSE).
