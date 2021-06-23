# number of images in each axis
Nim = [3, 3, 3]
# space group
space_group = 113

# charge, for Gaussian calculation
charge = 5
# spin multiplicity
sm = 1
# Gaussian route card
gau_memory = 120 # unit : GB
gau_nprocs = 28
gau_method = 'b3lyp'
gau_basis = '6-31g*'
gau_maxcycle = 4096
gau_conv = 6 # 8 is default
gau_symm = True

# convergence criteria of PCC
max_diff_crt = 1e-4

# print level
lcp2k = False
