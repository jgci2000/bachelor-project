import time
import numpy as np
from matplotlib import pyplot as plt

# Compute the therodynamics for the Ising Model
# João Inácio, 19th Jan. 2021

def main():
    start = time.process_time()
    
    # Some variables
    kB = 1
    Tc_Onsager = 2.269
    
    dim = "2D"
    lattice = "SS"
    NN = 4
    
    L = 8
    N_SPINS = 1 * L ** 2
    q_max = N_SPINS // 2 + 1
    REP = 10**6
    skip = N_SPINS
    
    run_max = 1000
    
    max_E = (1 / 2) * NN * N_SPINS
    max_M = N_SPINS

    NE = int(1 + (max_E / 2))
    NM = N_SPINS + 1
    NT = 100
    
    energies = np.linspace(- max_E, max_E, NE)
    magnetizations = np.linspace(- max_M, max_M, NM)
    
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
            
    # Energies
    E = np.zeros(len(temperatures))
    E2 = np.zeros(len(temperatures))
    mod_E = np.zeros(len(temperatures))
    
    for i in range(len(temperatures)):
        for j in range(len(energies)):
            E[i] += energies[j] * sum(JDOS[j, :]) * np.exp(- beta_vals[i] * energies[j]) / Z[i]
            E2[i] += (energies[j]**2) * sum(JDOS[j, :]) * np.exp(- beta_vals[i] * energies[j]) / Z[i]
            mod_E[i] += np.abs(energies[j]) * sum(JDOS[j, :]) * np.exp(- beta_vals[i] * energies[j]) / Z[i]
    
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
        
    # Mean magnetic susceptability, mean heat capacity and mean entropy
    mean_X = np.zeros(len(temperatures))
    mean_C = np.zeros(len(temperatures))
    mean_S = np.zeros(len(temperatures))
    
    for i in range(len(temperatures)):
        mean_X[i] = (M2[i] - M[i]**2) * beta_vals[i]
        mean_C[i] = (E2[i] - E[i]**2) * beta_vals[i]
        mean_S[i] = (E[i] - min_F[i]) / temperatures[i]
        
    # Heat capacity
    C = np.zeros(len(temperatures))
    
    F_sd = np.zeros(len(temperatures))
    h = np.abs(temperatures[1] - temperatures[2])
    for i in range(1, len(temperatures) - 1):
        F_sd[i] = (min_F[i - 1] - 2 * min_F[i] + min_F[i + 1]) / h**2
        
    for i in range(len(temperatures)):
        C[i] = - temperatures[i] * F_sd[i]
    
    # Tc apprxomation
    h = np.abs(temperatures[1] - temperatures[2])
    M_min_F_fd = np.zeros(len(temperatures))
    mod_M_fd = np.zeros(len(temperatures))
    
    for i in range(1, len(temperatures) - 1):
        M_min_F_fd[i] = (M_min_F[i + 1] - M_min_F[i - 1]) / (2 * h)
        mod_M_fd[i] = (mod_M[i + 1] - mod_M[i - 1]) / (2 * h)
    
    Tc_M_min_F = temperatures[np.where(M_min_F_fd == min(M_min_F_fd))[0][0]]
    Tc_mod_M = temperatures[np.where(mod_M_fd == min(mod_M_fd))[0][0]]
    Tc = Tc_M_min_F
    
    # Outputs
    print("Method computation time: {:.2f}s".format(run_time))
    print("Tc_M_min_F [L{:d}]: {:.3f}".format(L, Tc_M_min_F))
    print("Tc_mod_M [L{:d}]: {:.3f}".format(L, Tc_mod_M))

    plt.style.use('seaborn-whitegrid')
    
    # plt.figure(1)
    # plt.plot(temperatures, mod_M / N_SPINS, '.-b')

    # plt.xlabel("T")
    # plt.ylabel("<|M|>")
    # plt.title("<|M|> as a function of T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    # plt.figure(2)
    # plt.plot(temperatures, E / N_SPINS, '.-b')
    
    # plt.xlabel("T")
    # plt.ylabel("<E>")
    # plt.title("<E> as a function of T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    # plt.figure(3)
    # for i in range(0, len(temperatures)):
    #     plt.plot(magnetizations / N_SPINS, F[:, i] / N_SPINS, '-b', lw=1)
    #     plt.plot(M_min_F[i] / N_SPINS, min_F[i] / N_SPINS, '.b', ms=7.5)
    
    # plt.xlabel("M")
    # plt.ylabel("F")
    # plt.title("F as a function of M and T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    # plt.figure(4)
    # plt.plot(temperatures / Tc, M_min_F / N_SPINS, '.-b')
    
    # plt.xlabel("T/Tc")
    # plt.ylabel("M minF")
    # plt.title("Magnetization for F minina as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    # plt.figure(5)
    # plt.plot(temperatures / Tc, C / N_SPINS)
    
    # plt.xlabel("T/Tc")
    # plt.ylabel("C")
    # plt.title("Heat Capacity per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    # plt.figure(6)
    # plt.plot(temperatures / Tc, mean_C / N_SPINS)
    
    # plt.xlabel("T/Tc")
    # plt.ylabel("<C>")
    # plt.title("Mean Heat Capacity per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(temperatures / Tc, mod_M / N_SPINS, '.-b')
    # axs[0, 0].xlabel("T")
    # axs[0, 0].ylabel("<|M|>")
    axs[0, 0].set_title("<|M|> as a function of T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    axs[0, 1].plot(temperatures / Tc, E / N_SPINS, '.-b')
    # axs[0, 1].xlabel("T")
    # axs[0, 1].ylabel("<E>")
    axs[0, 1].set_title("<E> as a function of T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    axs[1, 0].plot(temperatures / Tc, M_min_F / N_SPINS, '.-b')
    # axs[1, 0].xlabel("T/Tc")
    # axs[1, 0].ylabel("M minF")
    axs[1, 0].set_title("Magnetization for F minina as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    for i in range(0, len(temperatures)):
        axs[1, 1].plot(magnetizations / N_SPINS, F[:, i] / N_SPINS, '-b', lw=1)
        axs[1, 1].plot(M_min_F[i] / N_SPINS, min_F[i] / N_SPINS, '.b', ms=7.5)
    # axs[1, 1].xlabel("M")
    # axs[1, 1].ylabel("F")
    axs[1, 1].set_title("F as a function of M and T | L = " + str(L) + " | REP = " + str(int(np.log10(REP))))
    
    fig, axs = plt.subplots(2, 2)
    
    axs[0, 0].plot(temperatures / Tc, C / N_SPINS, '.-b')
    # axs[0, 0].xlabel("T/Tc")
    # axs[0, 0].ylabel("C")
    axs[0, 0].set_title("Heat Capacity per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    axs[1, 0].plot(temperatures / Tc, mean_C / N_SPINS, '.-b')
    # axs[1, 0].xlabel("T/Tc")
    # axs[1, 0].ylabel("<C>")
    axs[1, 0].set_title("Mean Heat Capacity per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    axs[0, 1].plot(temperatures / Tc, mean_X / N_SPINS, '.-b')
    axs[0, 1].set_title("Mean Magnetic Susceptability per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    axs[1, 1].plot(temperatures / Tc, mean_S / N_SPINS, '.-b')
    axs[1, 1].set_title("Mean Entropy per spin as a function of T | L = " + str(L) + " | REP = 1E" + str(int(np.log10(REP))))
    
    print("Script runtime: {:.4f}s".format(time.process_time() - start))
    plt.show()
    
    
if __name__ == '__main__':
    main()
