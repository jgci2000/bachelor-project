#!/usr/bin/env python3

import sys

"""
int max_E = (1.0 / 2.0) * NN * N_atm;
int max_M = N_atm;
int NE = 1 + (max_E / 2);
int NM = N_atm + 1;
"""

L = int(sys.argv[1])
lattice = sys.argv[2].upper()
q_time_E = float(sys.argv[3])

if lattice == "SS":
	N_atm = L**2
	NN = 4
elif lattice == "SC":
	N_atm = L**3
	NN = 6
elif lattice == "BCC":
	N_atm = 2 * L**3
	NN = 8
elif lattice == "FCC":
	N_atm = 4 * L**3
	NN = 12
elif lattice == "HCP":
	N_atm = 2 * L**3
	NN = 12
elif lattice == "HEX":
	N_atm = L**3
	NN = 8

NE = 1 + 0.25 * NN * N_atm
NM = N_atm + 1

JDOS_dim = NE * NM
JDOS_nnz_half = JDOS_dim * 0.5 * 0.6

time_est = JDOS_nnz_half * q_time_E
print(f"Time estimate = {time_est}s; {time_est/60}min; {time_est/3600}h; {time_est/3600/24}days")



