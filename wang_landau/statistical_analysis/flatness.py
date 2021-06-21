#!/usr/bin/env python3

import numpy as np
from scipy.stats import norm
import sys

dim = "2D"
lattice = "SS"
NN = 4

L = 4
N_atm = 1 * L ** 2

max_E = (1 / 2) * NN * N_atm
max_M = N_atm
NE = int(1 + (max_E / 2))
NM = N_atm + 1
energies = np.linspace(- max_E, max_E, NE)
magnetizations = np.linspace(- max_M, max_M, NM)

q_max = (NM + 1) // 2 - 2
if NM % 2 == 0:
    q_max = NM // 2 - 3

flatness = int(sys.argv[1])

# Regular WL
f_final_exp_vals = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
f_final_vals = 1 + 10.0**(-f_final_exp_vals)

n_run = 100

JDOS_all = list()
for i in range(len(f_final_vals)):
    JDOS_all.append(list())

wall_time_all = np.zeros((n_run, len(f_final_vals)))

cfg_chkbrd = np.zeros((n_run, len(f_final_vals)))
cfg_slice = np.zeros((n_run, len(f_final_vals)))
cfg_zerozero = np.zeros((n_run, len(f_final_vals)))

JDOS_mean = list()
wall_time_mean = list()

for k, f_final_exp in enumerate(f_final_exp_vals):
    for run in range(1, n_run + 1):
        file_name = "".join(("./data/flatness/", str(flatness), "/", str(f_final_exp),"/",
                             str(run), "_JDOS_WL_Ising_", dim, "_", lattice, 
                             "_L", str(L), "_f", str(f_final_exp), 
                             "_flatness", str(flatness)))
        JDOS = np.loadtxt(file_name + ".txt")
        
        cfg_chkbrd[run - 1, k] = JDOS[len(energies) - 1, q_max+1]
        cfg_slice[run - 1, k] = JDOS[L, q_max+1]
        cfg_zerozero[run - 1, k] = JDOS[energies==0, magnetizations==0]

        JDOS_all[k].append(JDOS)

        with open(file_name + "_data.txt", 'r') as data_file:
            header = data_file.readline().strip("\n")
            
            for last_line in data_file:
                pass

            wall_time_all[run - 1, k] = float(last_line)
    
    JDOS_mean.append(sum(JDOS_all[k]) / n_run)
    wall_time_mean.append(sum(wall_time_all[:, k]) / n_run)

JDOS_exact = np.loadtxt('JDOS_exact_L4_SS.txt')
    
mean_error_abs = list()
mean_error = list()
var_error = list()

for k in range(len(f_final_exp_vals)):
    error_all = list()
    error_all_abs = list()
    
    for run in range(n_run):
        JDOS_error = JDOS_all[k][run] - JDOS_exact
        JDOS_error = JDOS_error[np.where(JDOS_exact > 0)[0], np.where(JDOS_exact > 0)[1]] / JDOS_exact[np.where(JDOS_exact > 0)[0], np.where(JDOS_exact > 0)[1]]
        error_all.append(np.sum(np.sum(JDOS_error)))

        JDOS_error_abs = np.abs(JDOS_all[k][run] - JDOS_exact)
        JDOS_error_abs = JDOS_error_abs[np.where(JDOS_exact > 0)[0], np.where(JDOS_exact > 0)[1]] / JDOS_exact[np.where(JDOS_exact > 0)[0], np.where(JDOS_exact > 0)[1]]
        error_all_abs.append(np.sum(np.sum(JDOS_error_abs)))
    
    mean_error_abs.append(np.mean(error_all_abs))
    
    fit_error = norm.fit(error_all)
    mean_error.append(fit_error[0])
    var_error.append(fit_error[1])

with open(f"mean_error_abs_L4_SS_p{flatness}.txt", 'w') as file:
    for i in range(len(f_final_vals)-2):
        file.write(f"{f_final_vals[i+2] - 1} {mean_error_abs[i+2]} {var_error[i+2]} \n")

