import time
import numpy as np
from matplotlib import pyplot as plot
from matplotlib import patches as mpatches

# Statistical Analisys script
# João Inácio, 19th Jan., 2021
# 
# TODO:
#   - runtime and time per E
#   

def main():
    start = time.process_time()
    
    kB = 1
    Tc_Onsager = 2.269
    
    dim = "2D"
    lattice = "SS"
    NN = 4
    run_max = 14
    
    L = 16
    N_SPINS = L * L
    q_max = N_SPINS // 2 + 1
    REP = 10**4
    skip = N_SPINS
    
    max_E = (1 / 2) * NN * N_SPINS
    max_M = N_SPINS

    NE = int(1 + (max_E / 2))
    NM = N_SPINS + 1
    
    energies = np.linspace(- max_E, max_E, NE)
    magnetizations = np.linspace(- max_M, max_M, NM)

    JDOS_avg = np.zeros((NE, NM))
    
    checkerboard_conf = np.zeros(run_max)
    slice_conf = np.zeros(run_max)
    zerozero_conf = np.zeros(run_max)

    for run in range(0, run_max):
        # JDOS = np.loadtxt("./Data/" + dim + "_" + lattice + "/JDOS_FSS_Ising_" + dim + "_" + lattice +
        #                 "_L" + str(L) + "_REP_1E" + str(int(np.log10(REP))) + "_skip_" + str(skip) + ".txt", usecols=range(NE))
        JDOS = np.loadtxt("16L/" + str(run + 1) + "_JDOS_2D_" + str(L) + "L_10E" + str(int(np.log10(REP))) + "_skip" + str(skip) + "_nRPS_CPP.txt")
        JDOS[:, q_max:NM] = JDOS[:, range(q_max-2, -1, -1)]
        
        checkerboard_conf[run] = JDOS[N_SPINS, N_SPINS // 2]
        for i in range(0, NE):
            if JDOS[i, N_SPINS // 2] != 0:
                slice_conf[run] = JDOS[i, N_SPINS // 2]
                break
        zerozero_conf[run] = JDOS[N_SPINS // 2, N_SPINS // 2]
        
        JDOS_avg += JDOS

    JDOS_avg /= run_max
    
    print(checkerboard_conf)
    print(slice_conf)
    print(zerozero_conf)
    

    print("Runtime:", time.process_time() - start)
    # plt.show()


if __name__ == '__main__':
    main()
    

