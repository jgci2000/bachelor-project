import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches

# Finite size effect for the Ising Model
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
    
    temperatures = np.linspace(0.1, 5, 1000, dtype=np.float128)
    betha_vals = 1 / (kB * temperatures)
    
    L_vals = np.array([4, 8, 16])
    Tc = np.zeros(len(L_vals))
    M2 = np.zeros((len(L_vals), len(temperatures)))
    M4 = np.zeros((len(L_vals), len(temperatures)))
    k = 0
    
    for L in L_vals:
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
        mod_M = np.zeros(len(temperatures))

        for i in range(0, len(temperatures)):        
            for q in range(0, NM):
                M[i] += magnetizations[q] * Z_M[q, i] / Z[i]
                M2[k, i] += (magnetizations[q]**2) * Z_M[q, i] / Z[i]
                M4[k, i] += (magnetizations[q]**4) * Z_M[q, i] / Z[i]
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
                Tc_L = temperatures[i]
                break
        
        Tc[k] = Tc_L
        k += 1
        print("Tc[{:d}]: {:f}".format(L, Tc_L))
    
    plt.figure(1)
    for i in range(0, len(Tc)):
        plt.plot(L_vals[i], Tc[i], 'ok')
    
    plt.plot([0, np.max(L_vals) + 6], [Tc_Onsager, Tc_Onsager], '-r')
    plt.xlim((0, np.max(L_vals) + 6))
    plt.ylim((Tc_Onsager - 0.5, np.max(Tc) + 0.5))
    
    plt.title("Tc vs L")
    black_pnts = mpatches.Patch(color='k', label='Tc')
    plt.legend(handles=[black_pnts])
    red_line = mpatches.Patch(color='r', label='Onsager Tc')
    plt.legend(handles=[red_line, black_pnts])
    
    plt.xlabel("L")
    plt.ylabel("Tc")
    
    plt.figure(2)
    inv_L_vals = 1 / L_vals
    for i in range(0, len(Tc)):
        plt.plot(inv_L_vals[i], Tc[i], 'ok')
    
    plt.plot([- 0.1, np.max(inv_L_vals) + 0.1], [Tc_Onsager, Tc_Onsager], '-r')
    plt.xlim((- 0.1, np.max(inv_L_vals) + 0.1))
    plt.ylim((Tc_Onsager - 0.5, np.max(Tc) + 0.5))
    
    a = np.polyfit(inv_L_vals, Tc, 1)
    
    x = np.linspace(0, 0.30, 100)
    y = a[0] * x + a[1]
    
    plt.plot(x, y, '-b')
    
    plt.title("Tc vs 1/L with linear regression (y={:.3f}x+{:.3f})".format(a[0], a[1]))
    black_pnts = mpatches.Patch(color='black', label='Tc')
    red_line = mpatches.Patch(color='r', label='Onsager Tc')
    blue_line = mpatches.Patch(color='b', label='LinReg')
    plt.legend(handles=[red_line, black_pnts, blue_line])
    
    plt.xlabel("1/L")
    plt.ylabel("Tc")
    
    L_Tc_Onsager = a[0] / (Tc_Onsager - a[1])
    print("L for Tc_Onsager ({:.3f}): {:d}".format(Tc_Onsager, int(L_Tc_Onsager)))
    
    Tc_Lin_Reg = a[1]
    print("Tc by the linear regression from fitted data ->", Tc_Lin_Reg)
    
    U = np.zeros((len(L_vals), len(temperatures)), dtype=np.float128)
    U = 1 - M4 / (3 * M2**2)
    plt.figure(3)
    for i in range(0, len(L_vals)):
        plt.plot(temperatures, U[i, :], '.-')

    plt.xlabel("T")
    plt.ylabel("U(T, L)")
    plt.title("Binber Cumulant as a function of L and T")
    
    idx = np.argwhere(np.diff(np.sign(U[1, :] - U[2, :]))).flatten()
    Tc_Binder = temperatures[idx[len(idx) - 1]]
    print("Tc for the Binder Cumulant ->", Tc_Binder)
    
    print("Runtime:", time.process_time() - start)
    plt.show()
    
    
if __name__ == '__main__':
    main()
