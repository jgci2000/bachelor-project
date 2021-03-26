import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches

# Finite size effect for the Ising Model
# João Inácio, 19th Jan. 2021

def main():
    start = time.process_time()
    
    # Some constants and global arrays
    kB = 1
    Tc_Onsager = 2.269
    
    dim = "2D"
    lattice = "SS"
    NN = 4

    NT = 100
    temperatures = np.linspace(0.1, 5, NT, dtype=np.float128)
    beta_vals = 1 / (kB * temperatures)
    
    L_vals = np.array([4, 8, 16])
    run_vals = np.array([1000, 1000, 50])
    REP_vals = np.array([10**6, 10**6, 10**4])
    
    # Initialization of thermodynamic arrays
    magnetizations = np.zeros(len(L_vals))
    M = np.zeros((len(L_vals), len(temperatures)))
    M2 = np.zeros((len(L_vals), len(temperatures)))
    M4 = np.zeros((len(L_vals), len(temperatures)))
    mod_M = np.zeros((len(L_vals), len(temperatures)))
    
    F = np.zeros((len(L_vals), len(temperatures)))
    M_min_F = np.zeros((len(L_vals), len(temperatures)))
    min_F = np.zeros((len(L_vals), len(temperatures)))
    
    E = np.zeros((len(L_vals), len(temperatures)))
    E2 = np.zeros((len(L_vals), len(temperatures)))
    mod_E = np.zeros((len(L_vals), len(temperatures)))

    mean_X = np.zeros((len(L_vals), len(temperatures)))
    mean_C = np.zeros((len(L_vals), len(temperatures)))
    mean_S = np.zeros((len(L_vals), len(temperatures)))
    C = np.zeros((len(L_vals), len(temperatures)))

    Tc_M_min_F = np.zeros(len(L_vals))
    Tc_mod_M = np.zeros(len(L_vals))
    
    run_time_vals = np.zeros(len(L_vals))
    JDOS_vals = np.zeros(len(L_vals))
    
    for k in range(len(L_vals)):
        L = L_vals[k]
        N_SPINS = L * L
        q_max = N_SPINS // 2 + 1
        REP = REP_vals[k]
        skip = N_SPINS // 2
        
        run_max = run_vals[k]
        
        max_E = (1 / 2) * NN * N_SPINS
        max_M = N_SPINS

        NE = int(1 + (max_E / 2))
        NM = N_SPINS + 1
        
        energies = np.linspace(- max_E, max_E, NE)
        magnetizations[k] = np.linspace(- max_M, max_M, NM)
        
        print("System -> " + dim + "_" + lattice + " | L" + str(L) + " | REP: " + str(int(np.log10(REP))) + " | NT: " + str(NT))
        
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
        JDOS_all[k] = JDOS
        run_time_vals[k] = run_time / run_max
        
        # Partition function and Helmholtz free energy
        Z = np.zeros(len(temperatures), dtype=np.float128)
        Z_M = np.zeros((NM, len(temperatures)), dtype=np.float128)
        
        for q in range(0, NM):
            hits = np.where(JDOS[:, q] != 0)[0]

            for i in range(0, len(hits)):
                Z_M[q, :] += JDOS[hits[i], q] * np.exp(- beta_vals * energies[hits[i]])
            
            Z += Z_M[q, :]
            F[k][q, :] = - kB * temperatures * np.log(Z_M[q, :])
        
        # Magnetizations
        for i in range(0, len(temperatures)):        
            for q in range(0, NM):
                M[k][i] += magnetizations[k][q] * Z_M[q, i] / Z[i]
                M2[k][i] += (magnetizations[k][q]**2) * Z_M[q, i] / Z[i]
                M4[k][i] += (magnetizations[k][q]**4) * Z_M[q, i] / Z[i]
                mod_M[k][i] += np.abs(magnetizations[k][q]) * Z_M[q, i] / Z[i]
                
        # Energies
        for i in range(len(temperatures)):
            for j in range(len(energies)):
                E[k][i] += energies[j] * sum(JDOS[j, :]) * np.exp(- beta_vals[i] * energies[j]) / Z[i]
                E2[k][i] += (energies[j]**2) * sum(JDOS[j, :]) * np.exp(- beta_vals[i] * energies[j]) / Z[i]
                mod_E[k][i] += np.abs(energies[j]) * sum(JDOS[j, :]) * np.exp(- beta_vals[i] * energies[j]) / Z[i]
        
        # Magnetization for F minima
        for i in range(0, len(temperatures)):
            min_F[k][i] = F[0, i]
            q_min = 0
            
            for q in range(0, len(magnetizations[k])):
                if F[q, i] < min_F[k][i]:
                    min_F[k][i] = F[q, i]
                    q_min = q
            
            M_min_F[k][i] = np.abs(magnetizations[k][q_min])
            
        # Mean magnetic susceptability, mean heat capacity and mean entropy
        for i in range(len(temperatures)):
            mean_X[k][i] = (M2[k][i] - M[k][i]**2) * beta_vals[i]
            mean_C[k][i] = (E2[k][i] - E[k][i]**2) * beta_vals[i]
            mean_S[k][i] = (E[k][i] - min_F[k][i]) / temperatures[i]
            
        # Heat capacity
        F_sd = np.zeros(len(temperatures))
        h = np.abs(temperatures[1] - temperatures[2])
        for i in range(1, len(temperatures) - 1):
            F_sd[i] = (min_F[k][i - 1] - 2 * min_F[k][i] + min_F[k][i + 1]) / h**2
            
        for i in range(len(temperatures)):
            C[k][i] = - temperatures[i] * F_sd[i]
        
        # Tc apprxomation
        h = np.abs(temperatures[1] - temperatures[2])
        M_min_F_fd = np.zeros(len(temperatures))
        mod_M_fd = np.zeros(len(temperatures))
        
        for i in range(1, len(temperatures) - 1):
            M_min_F_fd[i] = (M_min_F[i + 1] - M_min_F[i - 1]) / (2 * h)
            mod_M_fd[i] = (mod_M[i + 1] - mod_M[i - 1]) / (2 * h)
        
        Tc_M_min_F[k] = temperatures[np.where(M_min_F_fd == min(M_min_F_fd))[0][0]]
        Tc_mod_M[k] = temperatures[np.where(mod_M_fd == min(mod_M_fd))[0][0]]
        
        # Normalize computations
        magnetizations[k] /= N_SPINS
        mod_M[k] /= N_SPINS
        M2[k] /= N_SPINS
        M4[k] /= N_SPINS
        E[k] /= N_SPINS
        M_min_F[k] /= N_SPINS
        F[k] /= N_SPINS
        min_F[k] /= N_SPINS
        C[k] /= N_SPINS
        mean_C[k] /= N_SPINS
        mean_X[k] /= N_SPINS
        mean_S[k] /= N_SPINS

    # Plots
    print()
    for k in range(len(L_vals)):
        print("L{:d} -> Tc_M_min_F: {:.3f}; Tc_mod_M: {:.3f}".format(L_vals[k], Tc_M_min_F[k], Tc_mod_M[k]))

    plt.style.use('seaborn-whitegrid')
    
    fig, axs = plt.subplots(2, 2)
    for k in range(len(L_vals)):
        axs[0, 0].plot(temperatures, mod_M[k], '.-b')
        # axs[0, 0].xlabel("T")
        # axs[0, 0].ylabel("<|M|>")
        axs[0, 0].set_title("<|M|> as a function of T")
        
        axs[0, 1].plot(temperatures, E[k], '.-b')
        # axs[0, 1].xlabel("T")
        # axs[0, 1].ylabel("<E>")
        axs[0, 1].set_title("<E> as a function of T")
        
        axs[1, 0].plot(temperatures, M_min_F[k], '.-b')
        # axs[1, 0].xlabel("T/Tc")
        # axs[1, 0].ylabel("M minF")
        axs[1, 0].set_title("Magnetization for F minina as a function of T")
        
        for i in range(0, len(temperatures)):
            axs[1, 1].plot(magnetizations, F[k][:, i], '-b', lw=1)
            axs[1, 1].plot(M_min_F[i], min_F[i], '.b', ms=7.5)
        # axs[1, 1].xlabel("M")
        # axs[1, 1].ylabel("F")
        axs[1, 1].set_title("F as a function of M and T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    fig, axs = plt.subplots(2, 2)
    
    axs[0, 0].plot(temperatures, C, '.-b')
    # axs[0, 0].xlabel("T/Tc")
    # axs[0, 0].ylabel("C")
    axs[0, 0].set_title("Heat Capacity per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    axs[1, 0].plot(temperatures, mean_C, '.-b')
    # axs[1, 0].xlabel("T/Tc")
    # axs[1, 0].ylabel("<C>")
    axs[1, 0].set_title("Mean Heat Capacity per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    axs[0, 1].plot(temperatures, mean_X, '.-b')
    axs[0, 1].set_title("Mean Magnetic Susceptability per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    axs[1, 1].plot(temperatures, mean_S, '.-b')
    axs[1, 1].set_title("Mean Entropy per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    plt.figure(3)
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
    
    plt.figure(4)
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
    plt.figure(5)
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
