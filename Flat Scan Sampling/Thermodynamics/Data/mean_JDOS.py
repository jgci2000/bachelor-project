import time
import numpy as np


def main():
    start = time.process_time()
    
    kB = 1
    Tc_Onsager = 2.269
    
    dim = "3D"
    lattice = "SC"
    NN = 6
    
    run_max = 1000
    
    L = 4
    N_SPINS = 1 * L ** 3
    q_max = N_SPINS // 2 + 1
    REP = 10**5
    skip = N_SPINS
    
    max_E = (1 / 2) * NN * N_SPINS
    max_M = N_SPINS

    NE = int(1 + (max_E / 2))
    NM = N_SPINS + 1
    
    energies = np.linspace(- max_E, max_E, NE)
    magnetizations = np.linspace(- max_M, max_M, NM)
    
    JDOS_mean = np.zeros((NE, NM))
    run_time_mean = np.zeros(1)
       
    for run in range(1, run_max + 1):
        file_name = "./" + dim + "_" + lattice + "/L" + str(L) + "/" + str(int(np.log10(REP))) + "/" + str(run) + "_JDOS_FSS_Ising_" + dim + "_" + lattice + "_L" + str(L) + "_REP_1E" + str(int(np.log10(REP))) + "_skip_" + str(skip)
                          
        JDOS = np.loadtxt(file_name + ".txt")
        JDOS[:, q_max:NM] = JDOS[:, range(q_max-2, -1, -1)]

        JDOS_mean += JDOS
        
        with open(file_name + "_data.txt", 'r') as data_file:
            header = data_file.readline().strip("\n")
            for i in range(0, q_max - 1):
                line = data_file.readline().strip("\n").split(" ")
            
            run_time_mean += float(data_file.readline().strip("\n"))
    
    JDOS_mean /= run_max
    run_time_mean /= run_max
    
    np.savetxt("./" + dim + "_" + lattice + "/mean_JDOS_FSS_Ising_" + dim + "_" + lattice + "_L" + str(L) + 
               "_REP_1E" + str(int(np.log10(REP))) + "_skip_" + str(skip) + ".txt", JDOS_mean, fmt='%.5e')
    np.savetxt("./" + dim + "_" + lattice + "/mean_data_FSS_Ising_" + dim + "_" + lattice + "_L" + str(L) + 
               "_REP_1E" + str(int(np.log10(REP))) + "_skip_" + str(skip) + ".txt", run_time_mean, fmt='%.5e')

    print("Runtime:", time.process_time() - start)
    
    
    
if __name__ == '__main__':
    main()

