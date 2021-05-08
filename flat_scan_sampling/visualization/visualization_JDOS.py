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

M, E = np.meshgrid(magnetizations, energies)

q = 8
snap_max = 171

JDOS_complete = np.loadtxt("./data_JDOS/L4_SS/" + str(snap_max) + "_JDOS.txt")

fig = plt.figure(1, figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

JDOS_all = list()
log_JDOS_all = list()

i = 0
for snapshot in range(0, snap_max + 1):
    JDOS = np.loadtxt("./data_JDOS/L4_SS/" + str(snapshot) + "_JDOS.txt")

    JDOS_dif_0 = np.where(JDOS > 0)
    
    log_JDOS_all.append(np.zeros((NE, NM)))
    log_JDOS_all[i][JDOS_dif_0[0], JDOS_dif_0[1]] = np.log10(JDOS[JDOS_dif_0[0], JDOS_dif_0[1]])
    
    JDOS_all.append(JDOS)
    i += 1

cmap_show = True
def update(snapshot):
    global cmap_show
    
    ax.clear()
    
    # surf = ax.plot_surface(M, E, log_JDOS_all[snapshot], vmin=np.min(log_JDOS_all[i-1]), vmax=np.max(log_JDOS_all[i-1]),
    #                       cmap=cm.cividis, antialiased=False, rstride=1, cstride=1)
    
    E_h = E.flatten()
    M_h = M.flatten()
    JDOS_h = np.zeros_like(E_h)

    dE = 2.5 * np.ones_like(E_h)
    dM = 1.25 * np.ones_like(E_h)
    dJDOS = log_JDOS_all[snapshot].flatten()
    
    cmap = cm.get_cmap('jet')
    max_height = np.max(log_JDOS_all[snap_max])
    min_height = np.min(log_JDOS_all[snap_max])
    rgba = [cmap((k-min_height)/max_height) for k in dJDOS] 

    bar = ax.bar3d(M_h, E_h, JDOS_h, dM, dE, dJDOS, color=rgba)
    
    if (cmap_show):
        fig.colorbar(bar, shrink=0.5, aspect=10)
        cmap_show = False
        
    # ax.set_ylim([0, 1.1 * np.log10(np.max(JDOS_complete))])
    # ax.set_title(f"Construction of the JDOS by FSS")
    # ax.set_ylabel("E")
    # ax.set_xlabel("M")
    # ax.set_zlabel("log10(JDOS)")
    
    ax.view_init(90, -90)
    ax.set_axis_off()

    print(str(snapshot) + "/" + str(i-1), end='\r')

ani = FuncAnimation(fig, update, np.arange(0, len(log_JDOS_all), 1))
# ani.save("animation_FSS_" + dim + "_" + lattice + "_JDOS.gif", fps=30)#, writer='imagemagick', fps=30)



plt.show()


