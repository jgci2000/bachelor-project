import time
import numpy as np
from matplotlib import pyplot as plt

# Compute the therodynamics for the Ising Model
# João Inácio, 19th Jan. 2021
# 
# TODO:
#   - Cv and energy stuff


def print_matrix(matrix, NE, NM):
    for i in range(0, NE):
        for j in range(0, NM):
            print(matrix[i, j], end=" ")
        print()

def main():
    start = time.process_time()
    
    kB = 1
    Tc_Onsager = 2.269
    
    dim = "2D"
    lattice = "SS"
    NN = 4
    
    L = 4
    N_SPINS = L * L
    q_max = N_SPINS // 2 + 1
    REP = 10**4
    skip = N_SPINS // 2
    
    max_E = (1 / 2) * NN * N_SPINS
    max_M = N_SPINS

    NE = int(1 + (max_E / 2))
    NM = N_SPINS + 1
    
    energies = np.linspace(- max_E, max_E, NE)
    magnetizations = np.linspace(- max_M, max_M, NM)

    JDOS = np.loadtxt("./Data/" + dim + "_" + lattice + "/JDOS_FSS_Ising_" + dim + "_" + lattice +
                        "_L" + str(L) + "_REP_1E" + str(int(np.log10(REP))) + "_skip_" + str(skip) + ".txt", usecols=range(NE))
    JDOS[:, q_max:NM] = JDOS[:, range(q_max-2, -1, -1)]
    
    temperatures = np.linspace(0.1, 5, 1000, dtype=np.float128)
    betha_vals = 1 / (kB * temperatures)
    
    Z = np.zeros(len(temperatures), dtype=np.float128)
    Z_M = np.zeros((NM, len(temperatures)), dtype=np.float128)
    F = np.zeros((NM, len(temperatures)))
        
    for q in range(0, NM):
        hits = np.where(JDOS[:, q] != 0)[0]

        for i in range(0, len(hits)):
            Z_M[q, :] += JDOS[hits[i], q] * np.exp(- betha_vals * energies[hits[i]])
        
        Z += Z_M[q, :]
        F[q, :] = - kB * temperatures * np.log(Z_M[q, :])
    
    M = np.zeros(len(temperatures))
    M2 = np.zeros(len(temperatures))
    mod_M = np.zeros(len(temperatures))

    for i in range(0, len(temperatures)):        
        for q in range(0, NM):
            M[i] += magnetizations[q] * Z_M[q, i] / Z[i]
            M2[i] += (magnetizations[q]**2) * Z_M[q, i] / Z[i]
            mod_M[i] += np.abs(magnetizations[q]) * Z_M[q, i] / Z[i]
    
    M_min_F = np.zeros(len(temperatures))
    
    for i in range(0, len(temperatures)):
        min_M = F[0, i]
        for q in range(0, len(magnetizations)):
            if F[q, i] < min_M:
                min_M = F[q, i]
        
        q = np.where(F[:, i] == min_M)[0]
        M_min_F[i] = np.abs(magnetizations[q[0]])
        
    for i in range(0, len(temperatures)):
        if M_min_F[i] == 0:
            Tc = temperatures[i]
            break
    
    print("Tc[{:d}]: {:f}".format(L, Tc))

    plt.figure(1)
    plt.plot(temperatures, mod_M / N_SPINS, '.-')
    plt.xlabel("T")
    plt.ylabel("<|M|>")
    plt.title("Mean Absolute magnetization as a function of T")
    
    plt.figure(2)
    for i in range(0, len(temperatures)):
        plt.plot(magnetizations / N_SPINS, F[:, i] / N_SPINS)
    plt.xlabel("M")
    plt.ylabel("F")
    plt.title("Helmholtz Free Energy as a function of M and T")
    
    plt.figure(3)
    plt.plot(temperatures / Tc, M_min_F / N_SPINS, '.-')
    plt.xlabel("T/Tc")
    plt.ylabel("M minF")
    plt.title("Magnetization for F minina as a function of temperature")
    
    print("Runtime:", time.process_time() - start)
    plt.show()
    
    
if __name__ == '__main__':
    main()
