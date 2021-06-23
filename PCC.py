import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from PCC_para import *

class PCC():

    def __init__(self):

        print('Polarized Crystal Charge (PCC)\n')
        print('Program Author : Wenbin FAN (fanwenbin@shu.edu.cn)')
        print('Thank Xiaoqing ZHU for the help of mastering PCC! \n')

        self.sym_delta = 1e-4
        p3elements = ['H', 'He',
                      'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                      'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
        self.p3elements = [i.upper() for i in p3elements]

        return

    def read_para(self):

        self.Nim = Nim
        self.charge = charge
        self.sm = sm
        self.space_group = space_group
        self.gau_memory = gau_memory
        self.gau_nprocs = gau_nprocs
        self.gau_method = gau_method
        self.gau_basis = gau_basis
        self.gau_maxcycle = gau_maxcycle
        self.gau_conv = gau_conv
        if gau_symm:
            self.gau_symm = ''
        else:
            self.gau_symm = 'nosymm'
        self.max_diff_crt = max_diff_crt
        self.lcp2k = lcp2k
        
        return

    # https://stackoverflow.com/a/13849249/71522
    def angle_between(self, v1, v2):
        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def const2para(self):

        va = self.lattice_const[:, 0] # vector a
        vb = self.lattice_const[:, 1]
        vc = self.lattice_const[:, 2]

        a, b, c = [np.linalg.norm(x) for x in [va, vb, vc]]

        alpha = self.angle_between(vb, vc)
        beta = self.angle_between(va, vc)
        gamma = self.angle_between(va, vb)

        alpha, beta, gamma = [x * 180 / np.pi for x in [alpha, beta, gamma]]

        self.lattice_para = [a,b,c, alpha, beta, gamma]

        return

    def print_resp_option(self):

        f = open('resp.in', 'w')
        f.write('7\n18\n5\n1\n../sym_atom\n1\n')

        elements = list(set(self.ele_list))

        for element in elements:
            if element.upper() not in self.p3elements:
                f.write('\n')

        f.write('y\n0\n0\nq\n')
        f.close()

        return

    def read_sym_operator(self, pg):

        print('Reading symmetry operator `{}`...'.format(pg))
        path = r'space_group/SO_{}'.format(pg)
        try:
            f = open(path, 'r').readlines()
        except:
            print('[INFO] No symmetry operator file. ')
            return

        Nop = int(f[1])
        op = np.zeros([3,3,Nop])
        sf = np.zeros([3,Nop])

        for i in range(Nop):
            for j in range(3):
                line = f[5*i + 3+j].strip().split()
                assert len(line) == 3, '[ERR] Wrong operation matrix! line: `{}`'.format(line)
                op[:,j,i] = [float(x) for x in line]
            line = f[5*i + 3 + 3].split()
            if len(line) == 3:
                sf[:, i] = [float(x) for x in line]
            elif len(line) == 0:
                sf[:, i] = 0.0
            else:
                print('[ERR] Wrong operation shift! line: `{}`'.format(line))
                print(line)
                exit()

        self.Nop = Nop
        self.op = op
        self.sf = sf

        print('Point group number : {}'.format(pg))
        print('Number of symmetry operator : {}'.format(Nop))

        return

    def do_judge_sym(self):

        print('Judging symmetry atoms ...')
        q = np.dot(self.lattice_inv, self.q)
        q -= 0.5 # shift to origin
        Natom = self.Natom

        op_mat = self.op
        sf_mat = self.sf
        Nop = self.Nop

        judge_record = np.zeros(Natom)

        f = open('sym_atom', 'w')
        if lcp2k:
            k = open('sym_atom_CP2K', 'w')

        for i in range(Natom-1):
            if judge_record[i] > -50:
                judge_record[i] = -10
            sym_list = [i]

            for j in range(i+1, Natom):
                if judge_record[j] > -50 and judge_record[i] > -50: # not judged atom
                    for op in range(Nop):
                        q_old = q[:,i]
                        q_new = np.dot(q[:,j], op_mat[:, :, op])
                        q_new += sf_mat[:, op]

                        half = 0.5 - self.sym_delta
                        for c in range(3):
                            if q_old[c] < -half:
                                q_old[c] += 1
                            elif q_old[c] > 0.5:
                                q_old[c] -= 1
                            if q_new[c] < -half:
                                q_new[c] += 1
                            elif q_new[c] > 0.5:
                                q_new[c] -= 1

                        # print(i+1, j+1, q_old, q_new)
                        q_dif = q_new - q_old
                        if np.dot(q_dif, q_dif) < 1e-3:
                            # print(i, j, op)
                            sym_list.append(j)
                            judge_record[j] = -100
                            break

            if len(sym_list) > 1:
                for s in sym_list:
                    f.write('{}, '.format(s+1))
                f.write('\n')
                
                if lcp2k:
                    k.write('&CONSTRAINT\n  EQUAL_CHARGES\n  ATOM_LIST')
                    for s in sym_list:
                        k.write(' {}'.format(s+1))
                    k.write('\n&END\n')

        return

    def read_vasp(self, path='CONTCAR'):
        if not os.path.exists(path):
            if os.path.exists('CONTCAR'):
                path = 'CONTCAR'
            elif os.path.exists('POSCAR'):
                path = 'POSCAR'
            else:
                print('[ERR] file not exists! `{}`'.format(path))
                exit()
        f = open(path, 'r').readlines()

        factor = float(f[1])
        lattice_const = np.zeros([3,3])

        lattice_const[:, 0] = [float(x) for x in f[2].split()]
        lattice_const[:, 1] = [float(x) for x in f[3].split()]
        lattice_const[:, 2] = [float(x) for x in f[4].split()]

        lattice_const *= factor

        # read element list
        ele_name = f[5].split()
        ele_num = [int(x) for x in f[6].split()]
        assert len(ele_num) == len(ele_name)
        Nele = len(ele_name)

        ele_list = []
        for ele in range(Nele):
            for i in range(ele_num[ele]):
                ele_list.append(ele_name[ele])
        assert len(ele_list) == np.sum(ele_num)
        Natom = len(ele_list)

        q = np.zeros([3, Natom])
        for atom in range(Natom):
            q[:, atom] = [float(x) for x in f[8+atom].split()[:3]]

        # convert to Cartesian coordinate
        if f[7].strip()[0].lower() == 'd':
            q = np.dot(lattice_const, q)

        self.Natom = Natom
        self.ele_list = ele_list
        self.q = q
        self.lattice_const = lattice_const
        self.lattice_inv = np.linalg.inv(lattice_const)

        # convert lattice constant to lattice parameter (a, b, c, alpha, beta, gamma)
        self.const2para()

        print('Number of atoms : {}'.format(Natom))
        print('Element list : {}'.format(ele_name))
        print('Lattice parameter (Angstrom and degree) : ')
        print('   a={:8.3f}       b={:8.3f}      c={:8.3f}'.format(*self.lattice_para[:3]))
        print('   alpha={:8.3f}   beta={:8.3f}   gamma={:8.3f}'.format(*self.lattice_para[3:]))
        print('')

        return

    def write_init_gau(self):

        wkdir = '0'
        if not os.path.exists(wkdir):
            os.mkdir(wkdir)

        g = open(os.path.join(wkdir, '0.gjf'), 'w')
        g.write('%chk=0.chk\n')
        g.write('%nprocs={}\n'.format(self.gau_nprocs))
        g.write('%mem={}GB\n'.format(self.gau_memory))
        g.write('#p {}/genecp '.format(self.gau_method))
        g.write('{} \n'.format(self.gau_symm))
        g.write('scf(maxcycle={}, conver={})\n'.format(self.gau_maxcycle, self.gau_conv))
        g.write('\n')
        g.write('initial ab inito calculation for PCC\n\n')
        g.write('{} {}\n'.format(self.charge, self.sm))

        for atom in range(self.Natom):
            g.write('{}    {:18.12f}    {:18.12f}    {:18.12f}\n'.format(self.ele_list[atom], *self.q[:,atom]))
        g.write('\n')
        g.write('-C -N -O -H -F -Cl 0\n{}\n****\n'.format(self.gau_basis))
        g.write('-Mn -Br -I 0\nSDD\n****\n\n')
        g.write('-Mn -Br 0\nSDD\n\n')
        # g.write('antechamber-ini.esp\n\nantechamber.esp\n\n')
        g.write('\n\n\n')

        g.close()
        return

    def judge_sym(self):

        self.read_sym_operator(self.space_group)
        self.do_judge_sym()

        return

    def read_chg(self, it):

        print('Reading charge file in last iteration `{}`...'.format(it))
        chg_file = os.path.join(str(it), '{}.chg'.format(it))
        assert os.path.exists(chg_file), '[ERR] The requested charge file not exists! '

        f = open(chg_file, 'r').readlines()
        assert len(f) == self.Natom

        # chg = np.zeros(Natoms)
        chg = [float(line.split()[-1]) for line in f]
        assert len(chg) == self.Natom

        self.chg = chg

        return

    def get_pbc_list(self, im):

        l = []
        if im == 0:
            return [0]
        else:
            for i in range(im+1):
                if i != 0:
                    l.append(-i)
                    l.append(i)
                else:
                    l.append(0)
        l.sort()

        return l

    def super_charge(self, it):

        print('Supercell, generate next Gaussian input file ...')
        # number of images in each axis (x, y, z)
        Nim = self.Nim # A squared box is recommended, which reads a~b~c.
        Nx = self.get_pbc_list(Nim[0])
        Ny = self.get_pbc_list(Nim[1])
        Nz = self.get_pbc_list(Nim[2])

        g = open(os.path.join(str(it), '{}.gjf'.format(it)), 'w') # Gaussian input
        p = open(os.path.join(str(it), '{}.pdb'.format(it)), 'w') # .pdb for view charge in VMD

        # write file header
        g.write('%oldchk=../{}/{}.chk\n'.format(it-1,it-1))
        g.write('%chk={}.chk\n'.format(it))
        g.write('%nprocs={}\n'.format(self.gau_nprocs))
        g.write('%mem={}GB\n'.format(self.gau_memory))
        g.write('#p {}/genecp '.format(self.gau_method))
        g.write('{} '.format(self.gau_symm))
        g.write('charge\n')
        g.write('scf(maxcycle={}, conver={}) '.format(self.gau_maxcycle, self.gau_conv))
        g.write('guess=read\n')
        g.write('\n')
        g.write('Iteration {}, ab inito calculation for PCC\n\n'.format(it))
        g.write('{} {}\n'.format(self.charge, self.sm))

        p.write('TITLE      Preview for calculating Polarized Crystal Charge (PCC)\n')
        p.write('REMARK     Program author : Wenbin FAN (fanwenbin@shu.edu.cn)\n')
        p.write('REMARK     Iteration {}\n'.format(it))
        p.write('CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} {}\n'.format(*self.lattice_para, self.space_group))

        # write coordinates
        for atom in range(self.Natom):
            g.write('{}    {:18.12f}    {:18.12f}    {:18.12f}\n'.format(self.ele_list[atom], *self.q[:,atom]))
            p.write('{:6s}{:5d} {:^4s}'.format('ATOM', atom+1, self.ele_list[atom]))
            p.write('{:1s}{:3s} {:1s}{:4d}{:1s}   '.format('', 'MOL', '', 1, ''))
            p.write('{:8.3f}{:8.3f}{:8.3f}'.format(*self.q[:,atom]))
            p.write('{:6.2f}{:6.2f}          '.format(1.0, self.chg[atom]))
            p.write('{:>2s}{:2s}\n'.format(self.ele_list[atom], ''))
        g.write('\n')
        p.write('TER     {:>3d}      MOL\n'.format(0))

        pbc_count = 1
        atom_count = 1
        q_all = self.q

        for x in Nx:
            for y in Ny:
                for z in Nz:
                    # print(x,y,z)
                    if [x, y, z] != [0, 0, 0]:
                        shift_vec = np.dot(self.lattice_const, [x,y,z])
                        for atom in range(self.Natom):
                            q_new = self.q[:, atom] + shift_vec
                            #check collision
                            q_check = np.add(q_all, -q_new[:, None])
                            q_check = np.sum(np.abs(q_check), axis=0)

                            if not np.min(q_check) < 1e-4:
                                q_all = np.append(q_all, np.transpose([q_new]), axis=1)
                                g.write('{:18.12f}    {:18.12f}    {:18.12f}    {:16.10f}\n'.format(*q_new, self.chg[atom]))
                                p.write('{:6s}{:5d} {:^4s}'.format('ATOM', atom_count, self.ele_list[atom]))
                                atom_count += 1
                                # pbc_count*self.Natom+atom+1
                                p.write('{:1s}{:3s} {:1s}{:4d}{:1s}   '.format('', 'CHG', '', 1, ''))
                                p.write('{:8.3f}{:8.3f}{:8.3f}'.format(*q_new))
                                p.write('{:6.2f}{:6.2f} {:>4d}     '.format(1.0, self.chg[atom], atom+1))
                                p.write('{:>2s}{:2s}\n'.format('X', ''))
                            # else:
                            #     print(x,y,z,atom,self.ele_list[atom])
                        p.write('TER     {:>3d}      MOL\n'.format(pbc_count))
                        pbc_count += 1

        g.write('\n')
        g.write('-C -N -O -H -Cl 0\n{}\n****\n'.format(self.gau_basis))
        g.write('-Mn -Br 0\nSDD\n****\n\n')
        g.write('-Mn -Br 0\nSDD\n\n')
        g.write('\n\n\n')
        p.write('END')

        return

    def conv_judeg(self, it):

        self.read_chg(it-2)
        chg_old = np.array(self.chg)
        self.read_chg(it-1)
        chg_new = np.array(self.chg)
        assert len(chg_old) == len(chg_new), '[ERR] Different shape of charge file! Iteration : {}'.format(it)

        chg_diff = np.abs(chg_new - chg_old)
        diff_max = np.max(chg_diff)
        # diff_var = np.var(chg_diff) # not used

        if diff_max < self.max_diff_crt:
            print('PCC converged! ')
            exit()
        else:
            print('PCC not converged! ')
            print('Max difference : {:12.8f}'.format(diff_max))

        return

    def init_calc(self):

        self.write_init_gau()
        self.judge_sym()
        self.print_resp_option()

        return

    def iter_calc(self, it=1):

        if it > 1:
            self.conv_judeg(it)

        if not os.path.exists(str(it)):
            os.mkdir(str(it))

        self.read_chg(it-1)
        self.super_charge(it)

        self.print_resp_option()
        self.plot_chg()

        return

    def plot_chg(self):

        it_list = []
        for root, dirs, files in os.walk('.'):
            for dir in dirs:
                if dir.isnumeric():
                    chg_path = os.path.join(dir, '{}.chg'.format(dir))
                    if os.path.exists(chg_path):
                        it_list.append(dir)
        Nit = len(it_list)
        it_list.sort(key=int)
        it_int = [int(k) for k in it_list]

        chg = np.zeros([self.Natom, Nit])
        chg_dif = np.zeros([self.Natom, Nit-1])

        for i, it in enumerate(it_list):
            chg_path = os.path.join(it, '{}.chg'.format(it))
            chg_file = open(chg_path, 'r').readlines()
            chg[:, i] = [float(line.split()[-1]) for line in chg_file]
            if i > 0:
                chg_dif[:, i-1] = np.abs(chg[:, i] - chg[:, i-1])

        for i in range(self.Natom):
            plt.scatter(it_int, chg[i, :])
            plt.plot(it_int, chg[i, :])
        plt.savefig('charge_evolution.svg', format='svg')
        plt.close()

        for i in range(self.Natom):
            plt.plot(it_int[1:], chg_dif[i, :])
        plt.yscale('log')
        plt.savefig('charge_evolution_difference.svg', format='svg')

        return

    def main(self):

        self.read_para()
        self.read_vasp()

        opt = str(sys.argv[1]).strip()

        if opt == '0':
            self.init_calc()
        elif opt == 'p':
            self.plot_chg()
        else:
            self.iter_calc(int(opt))

        return

a = PCC()
a.main()