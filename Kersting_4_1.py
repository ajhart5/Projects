import numpy as np
import pdb

# This script computes the phase impedance matrix, sequence impedance matrix,
# and the neutral transformation matrix for problem 4.1 of 
# Distribution System Modeling, 4th Edition, Kersting
#
#
# Problem 4.1:
#
# Phase conductors: 556,500 26/7 ACSR
# Neutral conductor: 4/0 ACSR
# 
#    Size | Stranding | Material | DIAM (in.) | GMR (ft) | RES (Ohm/mile) | Capacity (A)
# 556,500 |      26/7 |     ACSR |      0.927 |   0.0313 |         0.1859 | 730
#     4/0 |       6/1 |     ACSR |      0.563 |  0.00814 |          0.592 | 340
#
#
# Example 4.1
#
# pc_res = 0.306  # Ohm/mile for phase conductor
# pc_gmr = 0.0244  # GMR for phase conductor in feet
# nc_res = 0.5920   # Ohm/mile for neutral conductor
# nc_gmr = 0.00814 # GMR for neutral conductor in feet
#
# d[1] = 29j        # Line 1 - Phase a
# d[2] = 2.5 + 29j  # Line 2 - Phase b
# d[3] = 7 + 29j    # Line 3 - Phase c
# d[4] = 4 + 25j    # Line 4 - neutral

def compute_self_impedance(r_i, GMR_i):

    # Modified Carson's Equation for self-impedance (Kersting, Eq 4.41)

    z_ii = r_i + 0.09530 + 0.12134*1j*(np.log(1/GMR_i) + 7.93402)

    return z_ii

def compute_mutual_impedance(D_ij):

    # Modified Carson's Equation for mutual impedance (Kersting, Eq 4.42)

    z_ij = 0.09530 + 0.12134*1j*(np.log(1/D_ij) + 7.93402)

    return z_ij

def format_complex_for_excel(z):
    real = np.round(np.real(z), 4)
    imag = np.round(np.imag(z), 4)
    if imag >= 0:
        return f"{real}+{imag}j"
    else:
        return f"{real}{imag}j"

# Define constants
# Practice problem 4.1
pc_res = 0.1859  # Ohm/mile for phase conductor
pc_gmr = 0.0313  # GMR for phase conductor in feet
nc_res = 0.592   # Ohm/mile for neutral conductor
nc_gmr = 0.00814 # GMR for neutral conductor in feet

# Specify each conductor position on the pole in Cartesian coordinates using
# complex number notation.
d = np.zeros(5, dtype=complex)

# Practice problem 4.1
d[1] = 29j        # Line 1 - Phase a
d[2] = 7 + 29j    # Line 2 - Phase b
d[3] = 2.5 + 29j  # Line 3 - Phase c
d[4] = 4 + 25j    # Line 4 - neutral

# Compute the distance between the positions
D = np.zeros((5,5), dtype=complex)
D[1,2] = D[2,1] = abs(d[1] - d[2])
D[2,3] = D[3,2] = abs(d[2] - d[3])
D[3,1] = D[1,3] = abs(d[3] - d[1])
D[1,4] = D[4,1] = abs(d[1] - d[4])
D[2,4] = D[4,2] = abs(d[2] - d[4])
D[3,4] = D[4,3] = abs(d[3] - d[4])

# The diagonal terms of the distance matrix are the GMRs of the phase and neutral conductors
D[1,1] = pc_gmr
D[2,2] = pc_gmr
D[3,3] = pc_gmr
D[4,4] = nc_gmr

# Compute the primitive impedance matrix
z = np.zeros((5,5), dtype=complex)
for i in range(1,5):
    for j in range(1,5):
        if i == j:
            if i == 4:  # neutral conductor
                z[i,j] = compute_self_impedance(nc_res, D[i,i])
            else:  # phase conductors
                z[i,j] = compute_self_impedance(pc_res, D[i,i])
        else:
            z[i,j] = compute_mutual_impedance(D[i,j])

z_ij = z[1:4,1:4]  # 3x3 submatrix for the phase conductors
z_in = z[1:4,4].reshape(3,1)    # 3x1 submatrix for the mutual impedance between phase and neutral
z_nn = z[4,4]      # self-impedance of the neutral conductor
z_nj = z[4,1:4]    # 1x3 submatrix for the mutual impedance between neutral and phase conductors

# Compute the phase impedance matrix using Kron reduction
z_abc = z_ij - z_in * (1/z_nn) * z_nj

# Compute the neutral transformation matrix
t_n = (-(1/z_nn) * z_nj).reshape(1,3)

# Compute the sequence impedance matrix
a_s = np.exp(1j*2*np.pi/3)  # operator for sequence transformation
A_s = np.array([[1, 1, 1],
                [1, a_s**2, a_s],
                [1, a_s, a_s**2]])  # sequence transformation matrix

z_012 = np.linalg.inv(A_s) @ z_abc @ A_s

# Print matrices in CSV format
print("z_abc matrix:")
for i in range(z_abc.shape[0]):
    row = [format_complex_for_excel(z_abc[i,j]) for j in range(z_abc.shape[1])]
    print(",".join(row))

print("\nz_012 matrix:")
for i in range(z_012.shape[0]):
    row = [format_complex_for_excel(z_012[i,j]) for j in range(z_012.shape[1])]
    print(",".join(row))

print("\nt_n matrix:")
for i in range(t_n.shape[0]):
    row = [format_complex_for_excel(t_n[i,j]) for j in range(t_n.shape[1])]
    print(",".join(row))