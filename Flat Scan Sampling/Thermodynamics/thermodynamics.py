import time
import numpy as np
from scipy.stats import norm
from scipy import interpolate
from matplotlib import pyplot as plt

# Compute the therodynamics for the Ising Model
# João Inácio, 19th Jan. 2021
# 
# TODO:
#   - Cv and energy stuff


def main():
    start = time.process_time()
    
    # Some variables
    kB = 1
    Tc_Onsager = 2.269
    
    dim = "2D"
    lattice = "SS"
    NN = 4
    
    L = 16
    N_SPINS = 1 * L ** 2
    q_max = N_SPINS // 2 + 1
    REP = 10**4
    skip = N_SPINS
    
    run_max = 50
    
    max_E = (1 / 2) * NN * N_SPINS
    max_M = N_SPINS

    NE = int(1 + (max_E / 2))
    NM = N_SPINS + 1
    NT = 50                         # Number of temperatures
    
    energies = np.linspace(- max_E, max_E, NE)
    magnetizations = np.linspace(- max_M, max_M, NM)
    
    # Read all JDOS files
    JDOS_mean = np.zeros((NE, NM))
    JDOS_all = list()
    run_time = 0
    
    for run in range(1, run_max + 1):
        file_name = "./Data/" + dim + "_" + lattice + "/L" + str(L) + "/" + str(int(np.log10(REP))) + "/" + str(run) + "_JDOS_FSS_Ising_" + dim + "_" + lattice + "_L" + str(L) + "_REP_1E" + str(int(np.log10(REP))) + "_skip_" + str(skip)
                          
        JDOS = np.loadtxt(file_name + ".txt")
        JDOS[:, q_max:NM] = JDOS[:, range(q_max-2, -1, -1)]

        JDOS_all.append(JDOS)
        
        with open(file_name + "_data.txt", 'r') as data_file:
            header = data_file.readline().strip("\n")
            for i in range(0, q_max - 1):
                line = data_file.readline().strip("\n").split(" ")
            
            run_time += float(data_file.readline().strip("\n"))
    
    JDOS = sum(JDOS_all) / run_max
    run_time /= run_max
    
    # Compute thermodynamics from mean JDOS
    temperatures = np.linspace(0.1, 5, NT, dtype=np.float128)
    beta_vals = 1 / (kB * temperatures)
    
    # Partition function and Helmholtz free energy
    Z = np.zeros(len(temperatures), dtype=np.float128)
    Z_M = np.zeros((NM, len(temperatures)), dtype=np.float128)
    F = np.zeros((NM, len(temperatures)))
    
    for q in range(0, NM):
        hits = np.where(JDOS[:, q] != 0)[0]

        for i in range(0, len(hits)):
            Z_M[q, :] += JDOS[hits[i], q] * np.exp(- beta_vals * energies[hits[i]])
        
        Z += Z_M[q, :]
        F[q, :] = - kB * temperatures * np.log(Z_M[q, :])
    
    # Magnetizations
    M = np.zeros(len(temperatures))
    M2 = np.zeros(len(temperatures))
    mod_M = np.zeros(len(temperatures))
    
    for i in range(0, len(temperatures)):        
        for q in range(0, NM):
            M[i] += magnetizations[q] * Z_M[q, i] / Z[i]
            M2[i] += (magnetizations[q]**2) * Z_M[q, i] / Z[i]
            mod_M[i] += np.abs(magnetizations[q]) * Z_M[q, i] / Z[i]
    
    # Magnetization for F minima
    M_min_F = np.zeros(len(temperatures))
    min_F = np.zeros(len(temperatures))
    
    for i in range(0, len(temperatures)):
        min_F[i] = F[0, i]
        q_min = 0
        
        for q in range(0, len(magnetizations)):
            if F[q, i] < min_F[i]:
                min_F[i] = F[q, i]
                q_min = q
        
        M_min_F[i] = np.abs(magnetizations[q_min])
    
    # Tc -> approximation
    # Second derivative = 0 (...)
    
    h = 0.0001
    temperatures_interp = np.arange(min(temperatures) + h, max(temperatures), h)
    M_min_F_interp = interpolate.interp1d(temperatures, M_min_F, kind='cubic')
    M_min_F_interp = M_min_F_interp(temperatures_interp)
    
    M_min_F_interp_sd = np.zeros(len(temperatures_interp))
    M_min_F_interp_fd = np.zeros(len(temperatures_interp))
       
    for i in range(1, len(temperatures_interp) - 1):
        M_min_F_interp_sd[i] = (M_min_F_interp[i + 1] - 2 * M_min_F_interp[i] + M_min_F_interp[i - 1]) / (h ** 2)
        M_min_F_interp_fd[i] = (M_min_F_interp[i + 1] - M_min_F_interp[i - 1]) / (2 * h)
    
    for i in range(1, len(temperatures_interp) - 1):
        if M_min_F_interp_fd[i - 1] > M_min_F_interp_fd[i] and M_min_F_interp_fd[i + 1] > M_min_F_interp_fd[i]:
            first_min = i
            break
    
    h = temperatures[2] - temperatures[1]
    M_min_F_fd = np.zeros(len(temperatures))
    for i in range(1, len(temperatures) - 1):
        M_min_F_fd[i] = (M_min_F[i + 1] - M_min_F[i - 1]) / (2 * h)
    
    print(M_min_F_interp_fd[i] / N_SPINS)
    print(temperatures_interp[i])
    print(min(M_min_F_interp_fd))
    
    # print(M_min_F_interp_sd)
    # plt.plot(temperatures_interp, M_min_F_interp_sd)
    # plt.figure(2)
    plt.plot(temperatures_interp, M_min_F_interp_fd / N_SPINS)
    plt.figure(2)
    plt.plot(temperatures, M_min_F_fd / N_SPINS)
    plt.show()
    # for i in range(1, len(temperatures) - 1):
    #     if M_min_F_sd[i - 1] != 0 and M_min_F_sd[i + 1] != 0 and M_min_F_sd[i] == 0:
    #         Tc = temperatures[i]
    #         break
    return
    # print(M_min_F)
    # for i in range(0, len(temperatures)):
    #     if M_min_F[i] == 0:
    #         Tc = temperatures[i]
    #         break
    
    # Get error bars for the M vs T and M_minF vs T
    if L == 4 and lattice == "SS":
        error_bar = np.zeros(len(temperatures))
        
        for i in range(0, len(temperatures)):
            M_min_F_error = np.zeros(run_max)
            
            for run in range(1, run_max + 1):
                JDOS = JDOS_all[run - 1]
                
                F_T = np.zeros(NM)
                Z_M_T = np.zeros(NM)
                Z_T = 0
                
                for q in range(0, NM):
                    hits = np.where(JDOS[:, q] != 0)[0]

                    for j in range(0, len(hits)):
                        Z_M_T[q] += JDOS[hits[j], q] * np.exp(- beta_vals[i] * energies[hits[j]])
                    
                    Z_T += Z_M_T[q]
                    F_T[q] = - kB * temperatures[i] * np.log(Z_M_T[q])
                
                min_M_T = F_T[0]
                for q in range(0, len(magnetizations)):
                    if F_T[q] < min_M_T:
                        min_M_T = F_T[q]
                
                q = np.where(F_T == min_M_T)[0]
                M_min_F_error[run - 1] = np.abs(magnetizations[q[0]])
            
            # print(M_min_F_error)
            (mu, sigma) = norm.fit(M_min_F_error)
            error_bar[i] = sigma
            # print(i)
    
    print("Tc[L{:d}]: {:.4f}".format(L, Tc))

    plt.style.use('seaborn-whitegrid')
    
    # plt.figure(1)
    # plt.plot(temperatures, mod_M / N_SPINS, '.-')

    # plt.xlabel("T")
    # plt.ylabel("<|M|>")
    # plt.title("Mean Absolute magnetization as a function of T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    plt.figure(2)
    for i in range(0, len(temperatures)):
        plt.plot(magnetizations / N_SPINS, F[:, i] / N_SPINS, '.-b', ms=5)
        plt.plot(M_min_F[i] / N_SPINS, min_F[i] / N_SPINS, '.b', ms=7.5)
    
    plt.xlabel("M")
    plt.ylabel("F")
    plt.title("Helmholtz Free Energy as a function of M and T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    plt.figure(3)
    if L != 4 or lattice != "SS":
        plt.plot(temperatures / Tc, M_min_F / N_SPINS, '.-b')
    else:
        plt.errorbar(temperatures / Tc, M_min_F / N_SPINS, yerr=error_bar, fmt='.-b')
    
    plt.xlabel("T/Tc")
    plt.ylabel("M minF")
    plt.title("Magnetization for F minina as a function of temperature | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    print("Runtime: {:.4f}".format(time.process_time() - start))
    plt.show()
    
    
if __name__ == '__main__':
    main()
