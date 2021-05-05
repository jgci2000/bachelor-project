"""
    JDOS contruction visualization script for FSS method
    João Inácio, May 2nd, 2021
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation, PillowWriter  

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

q = 8
snap_max = 632

JDOS_complete = np.loadtxt("./data/L4_SS/" + str(q) + "_" + str(snap_max) + "_JDOS.txt")
hist_complete = np.loadtxt("./data/L4_SS/" + str(q) + "_" + str(snap_max) + "_hist.txt")

fig = plt.figure(1, figsize=(10, 10))
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

JDOS_all = list()
hist_all = list()

for snapshot in range(snap_max + 1):
    JDOS = np.loadtxt("./data/L4_SS/" + str(q) + "_" + str(snapshot) + "_JDOS.txt")
    hist = np.loadtxt("./data/L4_SS/" + str(q) + "_" + str(snapshot) + "_hist.txt")
    JDOS_all.append(JDOS)
    hist_all.append(hist)

def update(snapshot):
    ax1.clear()
    ax2.clear()
    
    ax1.plot(energies, JDOS_all[snapshot], '-or')
    ax2.plot(energies, hist_all[snapshot], '-or')
    
    ax1.set_ylim([0, 1.1 * np.max(JDOS_complete)])
    ax1.set_title(f"Construction of the DOS at {q+1}")
    ax1.set_ylabel(f"DOS at {q+1}")
    ax1.set_xlabel("E")

    ax2.set_ylim([0, 1.1 * np.max(hist_complete)])
    ax2.set_title(f"Histogram of the accepted configurations at {q}")
    ax2.set_ylabel(f"Histogram at {q}")
    ax2.set_xlabel("E")
    
    print(str(snapshot) + "/" + str(snap_max), end='\r')

ani = FuncAnimation(fig, update, np.arange(0, len(JDOS_all), 1))
print("animation done!")
ani.save("animation_FSS_" + dim + "_" + lattice + "_DOS_" + str(q+1) + ".gif", fps=30)#, writer='imagemagick', fps=30)
print("animation saved!")

# plt.show()


