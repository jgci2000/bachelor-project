"""
    JDOS contruction visualization script for FSS method
    João Inácio, May 2nd, 2021
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation, PillowWriter 
from mpl_toolkits.mplot3d import Axes3D 

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

snap_max = 1380

JDOS_complete = np.loadtxt("./data/L4_SS/" + str(snap_max) + "_JDOS.txt")

fig = plt.figure(1, figsize=(10, 10))
ax1 = fig.add_subplot(projection='3d')

JDOS_all = list()

for snapshot in range(snap_max + 1):
    JDOS = np.loadtxt("./data/L4_SS/" + str(snapshot) + "_JDOS.txt")
    JDOS_all.append(JDOS)

M, E = np.meshgrid(magnetizations, energies)

def update(snapshot):
    ax1.clear()
    
    ax1.bar3d(M, E, JDOS_all[snapshot], 1, 1)
    
    ax1.set_ylim([0, 1.1 * np.max(JDOS_complete)])
    ax1.set_title(f"Construction of the DOS at {1}")
    ax1.set_ylabel(f"DOS at {1}")
    ax1.set_xlabel("E")
    
    print(str(snapshot) + "/" + str(snap_max), end='\r')

ani = FuncAnimation(fig, update, np.arange(0, len(JDOS_all), 1))
print("animation done!")
# ani.save("animation_FSS_" + dim + "_" + lattice + "_DOS_" + str(1) + ".gif", fps=30)#, writer='imagemagick', fps=30)
print("animation saved!")

plt.show()


