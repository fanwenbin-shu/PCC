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

# Procedure

The python program and shell script will do these things to obtain a converged PCCs. 

1. Read `CONTCAR` or `POSCAR` in current folder `./`. Create folder `./0/`. Convert to Gaussian input file (`./0/0.gjf`). 
2. Find all symmetry atoms based on SO. All symmetry atoms will be written in one line at `./sym_atom`. This symmetry will be used in later RESP fitting. 
3. Gaussian SCF. The atomic charge will be fitted by Multiwfn. (Note that there are extra line break if there are heavy atoms. Because there are not vdW radius defined in Multiwfn. We will be prompted to confirm the default radius. Please check this if you are about to fit a big system. )
4. `<i>` represents the current step. Create folder `./<i>`. The point charge will be added. Gaussian input file is written in `./<i>/<i>.gjf`. The last `.chk` file will be read as initial orbit guess in order to accelerate the SCF. 
5. RESP fit. This procedure terminates when the maximum of charge difference between two steps is less than the value defined in `PCC_para.py` (usually `1E-6`). Or go back to 4th step and `i++`. 

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