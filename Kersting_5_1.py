import numpy as np
import pdb

# 5.1 Determine the phase admittance matrix and sequence admittance matrix 
# in μS/mile for the three-phase overhead line of Problem 4.1.
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

def compute_self_potential_coefficient(S_ii, RD_i):

    # Modified Carson's Equation for self-potential coefficient (Kersting, Eq 5.9)
    # S_ii is the distance from the conductor i to its image i' (ft)
    # RD_i is the radius of the conductor i (ft)
    # Output is in units of mile/microF
    
    P_ii = 11.17689*np.log(S_ii/RD_i)

    return P_ii

def compute_mutual_potential_coefficient(S_ij, D_ij):

    # Modified Carson's Equation for mutual potential coefficient (Kersting, Eq 5.10)
    # S_ij is the distance from the conductor i to the image of conductor j (ft)
    # D_ij is the distance between conductor i and conductor j (ft)
    # Output is in units of mile/microF
    
    P_ij = 11.17689*np.log(S_ij/D_ij)

    return P_ij

def format_complex_for_excel(z):
    real = np.round(np.real(z), 4)
    imag = np.round(np.imag(z), 4)
    if imag >= 0:
        return f"{real}+{imag}j"
    else:
        return f"{real}{imag}j"

##########################################################################################

# Specify each conductor position on the pole in Cartesian coordinates using
# complex number notation.
d = np.zeros(5, dtype=complex)

# # Practice problem 4.1
# d[1] = 29j        # Line 1 - Phase a
# d[2] = 7 + 29j    # Line 2 - Phase b
# d[3] = 2.5 + 29j  # Line 3 - Phase c
# d[4] = 4 + 25j    # Line 4 - neutral

# pc_res = 0.1859  # Ohm/mile for phase conductor
# pc_gmr = 0.0313  # GMR for phase conductor in feet
# pc_radius = 0.927/2/12  # radius of phase conductor in feet
# nc_res = 0.592   # Ohm/mile for neutral conductor
# nc_gmr = 0.00814 # GMR for neutral conductor in feet
# nc_radius = 0.563/2/12  # radius of neutral conductor in feet

# Example 4.1
d[1] = 29j        # Line 1 - Phase a
d[2] = 2.5 + 29j  # Line 2 - Phase b
d[3] = 7 + 29j    # Line 3 - Phase c
d[4] = 4 + 25j    # Line 4 - neutral

pc_res = 0.306  # Ohm/mile for phase conductor
pc_gmr = 0.0244  # GMR for phase conductor in feet
pc_radius = 0.721/2/12  # radius of phase conductor in feet
nc_res = 0.5920   # Ohm/mile for neutral conductor
nc_gmr = 0.00814 # GMR for neutral conductor in feet
nc_radius = 0.563/2/12  # radius of neutral conductor in feet

# Compute the distance between the positions
D = np.zeros((5,5), dtype=complex)
D[1,2] = D[2,1] = abs(d[1] - d[2])
D[2,3] = D[3,2] = abs(d[2] - d[3])
D[3,1] = D[1,3] = abs(d[3] - d[1])
D[1,4] = D[4,1] = abs(d[1] - d[4])
D[2,4] = D[4,2] = abs(d[2] - d[4])
D[3,4] = D[4,3] = abs(d[3] - d[4])
D[1,1] = pc_gmr # diagonal terms of the matrix are the GMRs of the conductors
D[2,2] = pc_gmr
D[3,3] = pc_gmr
D[4,4] = nc_gmr

# Compute the image distance matrix
S = np.zeros((5,5))
S[1,2] = S[2,1] = abs(d[1] - np.conj(d[2]))
S[2,3] = S[3,2] = abs(d[2] - np.conj(d[3]))
S[3,1] = S[1,3] = abs(d[3] - np.conj(d[1]))
S[1,4] = S[4,1] = abs(d[1] - np.conj(d[4]))
S[2,4] = S[4,2] = abs(d[2] - np.conj(d[4]))
S[3,4] = S[4,3] = abs(d[3] - np.conj(d[4]))
S[1,1] = abs(d[1] - np.conj(d[1]))
S[2,2] = abs(d[2] - np.conj(d[2]))
S[3,3] = abs(d[3] - np.conj(d[3]))
S[4,4] = abs(d[4] - np.conj(d[4]))

# Compute the primitive potential coefficient matrix
P_prim = np.zeros((5,5))
for i in range(1,5):
    for j in range(1,5):
        if i == j: # self
            if i == 4:  # neutral conductor
                P_prim[i,j] = compute_self_potential_coefficient(S[i,j], nc_radius)
            else:  # phase conductors
                P_prim[i,j] = compute_self_potential_coefficient(S[i,j], pc_radius)
        else: # mutual
            P_prim[i,j] = compute_mutual_potential_coefficient(S[i,j], D[i,j])

# Compute the reduced potential coefficient matrix
P_ij = P_prim[1:4,1:4]  # 3x3 submatrix for the phase conductors
P_in = P_prim[1:4,4].reshape(3,1)    # mutual potential coefficient between phase and neutral
P_jn = P_prim[4,1:4]    # mutual potential coefficient between neutral and phase conductors
P_nn = P_prim[4,4]      # self-potential coefficient of the neutral conductor

P_abc = P_ij - P_in * P_jn / P_nn

# Compute the capacitance matrix
C_abc = np.linalg.inv(P_abc)

# Compute the three-phase shunt admittance matrix
y_abc = 376.9911j * C_abc  # in μS/mile

# Compute the three phase sequence admittance matrix
a_s = np.exp(1j*2*np.pi/3)  # operator for sequence transformation
A_s = np.array([[1, 1, 1],
                [1, a_s**2, a_s],
                [1, a_s, a_s**2]])  # sequence transformation matrix

y_012 = np.linalg.inv(A_s) @ y_abc @ A_s # in μS/mile

# Print matrices in CSV format
print("y_abc matrix:")
for i in range(y_abc.shape[0]):
    row = [format_complex_for_excel(y_abc[i,j]) for j in range(y_abc.shape[1])]
    print(",".join(row))

print("\ny_012 matrix:")
for i in range(y_012.shape[0]):
    row = [format_complex_for_excel(y_012[i,j]) for j in range(y_012.shape[1])]
    print(",".join(row))