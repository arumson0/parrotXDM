# parrotXDM 🦜
This is a python tool to compute XDM dispersion energies, given the output of an electronic structure method already containing XDM data (e.g. Quantum ESPRESSO or FHI-aims). The use is to run with alternative damping functions and the Axilrod-Teller-Muto C9 dispersion term.

# USAGE:
The compiler adds the parrotxdm executable to your path. You may want to copy this to your .bashrc (Sorry windows users, you're on your own.)

`parrotxdm input.yaml`

The .yaml file follows this template:
```
input_path:
  - ./examples/ca2n_bulk_b86bpbe-xdmbj.scf.out
file_type: qe

functional: b86bpbe
damping: z
output_unit: Ry

run_settings:
  verbose_conv: true
  pairwise: true
  triples: true
  extrapolate_triples: true
  r2tol: 0.99
```
Tags:
`input_path`: path to the electronic structure method output containing the XDM data. Can be a list of multiple files.

`file_type`: specify whether the method is reading a `qe` or `aims` file.

`functional`: determines the damping parameters for XDM. Currently supports B86bPBE, PBE, and revPBE.

`damping`: type of damping function. Options are 'bj' or 'z'

`output_unit`: What unit you'd like the XDM energy to be expressed in. Options are Ha, Ry, eV, kcal/mol, kj/mol.

`override_a1`, `override_a2`: Options to manually specify a1 and a2 for BJ damping.

`override_zdamp`: Option to manually specify zdamp for Z-damping.

`run_settings`:

`verbose_conv`: prints the energy convergence over supercells if true.

`pairwise`: compute the pairwise dispersion contribution (C6, C8, C10)

`triples`: compute the ATM dispersion contribution

`c11`: (EXPERIMENTAL) compute the DDQ/DQD/QDD C11 dispersion energy contribution

`extrapolate_triples`: If true, when the changes in energy become log-linear with R^2 = r2tol, extrapolate the converged energy.

`r2tol`: defines the tolerance for the log-linear regression.

# How parrotXDM Works
Pyxdm is controlled by `driver.py`, which loads and calls all the python functions contained in parrotXDM. First, the driver reads the input yaml file with the `parse_config()` function in `/Input_Output/parse_yaml_input.py`. That is then loaded into a configuration for the run with the `Config` function in `/Input_Output/config.py`. This is where selections are set, such as the type of input file (QE or FHI-aims), choice of atomic pairs and/or triples, and the damping parameters. Default damping parameters are set based on the functional, and are read from `/Data/set_params.py`. Once the configuration for the run is set, parrotXDM reads the atomic positions and multipole moments from the DFT output file(s) specified in the input yaml. The functions to read these properties are found in `/Input_Output/readers.py`, which in turn read atomic data (e.g. reference polarizabilities and atomic numbers) from `/Data/atom_info.py`. Because parrotXDM was originally written with QE in mind, atomic positions are stored in terms of alat, rather than as fractional coordinates.

Computing Cn and the resulting energy contribution for pairs and triples follows the same workflow. First, the Cn coefficients are computed from the multipole moments and polarizabilities. This is done with functions in the files `Physics/cpair.py` or `Physics/ctriple.py`. The resulting Cn are printed by the driver. For pairwise, the critical radius and Zinv = 1/(Zi+Zj) (for Becke-Johnson and Z-damping, respectively) are also printed. These quantities are needed for both pairs and triples, so `cpair()` will always be called. Next, the energy is computed from `/Physics/epair.py` or `Physics/etriple.py`. Note that the leading character (e or c) for files in `/Physics/` helps to designate the file as being for coefficients or energy. The `epair()` function is responsible for converging the energy with respect to supercells, but it contains a call to a function called `cpair_kernel.pairwise_xdm()`. This function is a compiled Fortran subroutine, which computes the dispersion energy between atomic pairs. Fortran is used because of the nonuply-nested for-loops required for periodic three-body calculations; python would be too slow. The same infrastructure is used for pairs for consistency. The energy subroutine is called iteraively, increasing the supercell shell size until the energy difference converges to within 3.675E-8 Hartree, adding the energy contribition from each shell to a variable called `exdm_total`, which is returned by the `epair()` function. If `verbose_conv` is requested in the yaml, the energy convergence over supercells will be printed. The `etriple()` function in `/Physics/etriple.py` works in the same way as `epair()`, but the scaling is much worse due to the nonuply-nested for-loop. For this reason, parrotXDM can extrapolate the ATM contribution over supercells as log-linear. See the "Notes on the log-linear regression method" section of this document for details on accuracy. This extrapolation is handled in the `etriple()` function. `driver.py` prints the pairwise and triplewise XDM energies separately, along with the total XDM energy. If only pairs or triple are requested, then the total will match that quantity.

Currently, parrotXDM can optionally compute the C11 coefficients and their contribution to the disperison energy. This feature is HIGHLY experimental, will require its own parameterization, which has not been done yet. It should generally not be used.

# Notes on the log-linear regression method
The total energy for the ATM term can be extrapolated using a log-linear regression. Doing so can dramatically decrease the computational time required to compute the ATM contribution. Seven points with R^2 >= r2tol are required. This feature is experimental, but a small benchmark on nine simple crystalline system (mix of atomic, molecular, and metallic) suggests the following performance trend.
```
ATM DISPSERION ENERGY
R^2    MAE       ME
0.900: 0.003667% -0.003664%
0.950: 0.003667% -0.003664%
0.975: 0.003351% -0.003349%
0.990: 0.002897% -0.002881%
0.995: 0.000323%  0.000043%
0.999: 0.000323%  0.000043%
N = 9
```
The very small errors suggest that the method is generally robust. Based on this benchmark, r2tol is reccomended to be between 0.99 and 0.995.
