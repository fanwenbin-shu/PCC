# ePCC

A general python code for evaluating **polarized crystal charge** (PCC).

Author : **Wenbin FAN** (fanwenbin@shu.edu.cn)
Supervisor : Prof. **Yongle LI** (yongleli@shu.edu.cn)

The author thank Xiaoqing ZHU (carolin118@shu.edu.cn) for the understanding of the principle of PCC. 

Please CITE the article (10.1016/j.cplett.2019.06.019) if PCC or this ePCC code was used in your work. 

# Motive

The short-ranged electric interaction is significant in a MD simulation. Nearly all types of atomic charge only describe the charge of molecules in vacuum. But how does the atomic charge changed in the periodic crystal? 

Prof. Yongle proposed a refined method, polarized crystal charge (PCC), based on polarized protein-specific charge (PCC). The PCC is suitable for any periodic system. We can obtain PCCs easily using Gaussian and Multiwfn. 

# Usage

The Gaussian (>= 16, https://gaussian.com) and Multiwfn (>= 3.7, http://sobereva.com/multiwfn) is required. 

1. Prepare the P1 structure file `CONTCAR` in VASP format. The free visualization software VESTA (http://jp-minerals.org/vesta/en/) is recommend. 
2. Set the parameters in `PCC_para.py`. 
3. Define the space group and corresponding symmetry operation (SO) in `./space_group/SO_<space group number>`. 
4. Run `./batch.sh`, then the PCC calculation will start. 

## Parameters

Please modify the parameter file `PCC_para.py` before production run. The order of parameters does not effects. 

1. `Nim`, a one-dimensional integer list with length 3, means the images in each axis. `Nim = [3, 3, 3]` means that there will be $(2N_{\mathrm{im}}(0)+1)(2N_{\mathrm{im}}(1)+1)(2N_{\mathrm{im}}(2)+1) = 7 \times 7 \times 7 = 343$ images. 
2. `space_group`, integer. The SO file in `./space_group/SO_<space_group>`. Format is listed below. 
3. DFT calculation parameters. 
   1. `charge` and `sm` are two integers in Gaussian input file (.gjf), corresponding net charge and spin multiplicity of current system. Note that the net charge could not be zero probably if `POSCAR` include extra molecules. 
   2. `gau_memory` and `gau_nprocs` are memory amount (GB, integer) and number of processors (integer) in Gaussian DFT calculations. They will be written in the head of .gjf file as `%mem=<gau_memory>GB` and `%nprocs=<gau_nprocs>`. 
   3. `gau_method` and `gau_basis` are theoretical method and basis set (both characters) for Gaussian DFT calculation. They method will be written in the keyword line of .gjf file as `<gau_method>/genecp`. And the basis set will be written in the end of the .gjf file as `-C -H -O -N 0\n<gau_basis>`. Some heavy atoms (currently Mn and Br) use SDD pseudopotential and basis sets. 
   4. `gau_maxcycle` and `gau_conv` set the SCF parameter, maximum cycles and convergence criteria (both integers). They will be written in the keyword line of .gjf file as `SCF(maxcycle=<gau_maxcycle>, conver=<gau_conv>)`. Convergence criteria is 8 in default but it is too strict for SCF-only calculation especially for big system. So `gau_conv=6` is recommended. Please see 2(6) in "The way to solve non-convergence for DFT" (Chinese, http://sobereva.com/61) for more information. 
   5. `gau_symm` is a Boolean for symmetry. `True` for nothing and `False` for `nosymm` in keyword line. 
4. `max_diff_crt` is the criteria for PCC convergence (float). If the maximum difference of between two RESP charges for all atoms is less than `max_diff_crt`, the PCC calculation will terminate. `1e-5` is recommended for `gau_conv=8` and small system, and `1e-4` for `gau_conv=6` and bigger system. 
5. `lcp2k` is a old-fashioned option (Boolean). `True` for printing constraint in CP2K format. We no longer use CP2K in PCC calculation since the wave function is not compatible with Multiwfn (or need complicated manual processing) and the function of RESP charge is too simple to utilize. 

## Prepare initial configuration

For the non-periodic DFT calculation, a molecule must be held together. Some molecule in the edge of box should be repeated in each edge. For better visualization, the VESTA is recommended to join each part in the lattice cell. 

1. Load any supported format in VESTA. The structure with symmetry is recommended, like `.cif`. 
2. Delete extra molecules or add molecules by expanding the cell (supercell). All molecules must be held together. The total atoms must not be less than the atoms in original cell, otherwise the cell is not integrity. 
3. Export current structure in P1 format. All displayed atoms should be exported. The P1 format is the structure for VASP actually. 

## File format for symmetry operation

1. The first line is comment. Other lines will be read. 
2. The second line defines the number of SO, except the trivial operation (*x, y, z*). Blank line. 
3. Five lines are one SO. Three are symmetry matrix. One is translation vector. One blank line. 

Here is an example. The last SO of space group 113 is $(-y+\frac12,-x+\frac12,z)$, which equals in fractional coordinate
$$
\left(\matrix{x'\\y'\\z'}\right) = \left(\matrix{0&-1&0\\-1&0&0\\0&0&1}\right) \left(\matrix{x\\y\\z}\right) + \frac12\left(\matrix{1\\1\\0}\right).
$$
We should define the SO as below. 

```
 0 -1  0
-1  0  0
 0  0  1
0.5    0.5    0.0

```

# Program Procedure

The python program and shell script will do these things together to obtain a converged PCCs. 

1. Read `CONTCAR` or `POSCAR` in current folder `./`. Create folder `./0/`. Convert to Gaussian input file (`./0/0.gjf`). 
2. Find all symmetry atoms based on SO. All symmetry atoms will be written in one line at `./sym_atom`. This symmetry will be used in later RESP fitting. 
3. Gaussian SCF. The atomic charge will be fitted by Multiwfn. 
4. `<i>` represents the current step. Create folder `./<i>`. The point charge will be added. Gaussian input file is written in `./<i>/<i>.gjf`. The last `.chk` file will be read as initial orbit guess in order to accelerate the SCF. 
5. RESP fit. This procedure terminates when the maximum of charge difference between two steps is less than the value defined in `PCC_para.py` (usually `1E-6`). Or go back to 4th step and `i++`. 
