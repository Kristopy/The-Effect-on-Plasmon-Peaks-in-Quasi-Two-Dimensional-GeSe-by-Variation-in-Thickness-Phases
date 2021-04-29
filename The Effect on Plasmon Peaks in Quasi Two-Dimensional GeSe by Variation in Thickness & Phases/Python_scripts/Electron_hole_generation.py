import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})


T = np.linspace(0,273,1000)

def Generation_el_ho(T, E_g):
    k_b = 8.617333262145e-5
    q = 4.80320425e-10

    g = np.exp(-(q*E_g)/k_b*T)

    return g


g = Generation_el_ho(T, 0.2)

plt.figure()

plt.plot(T, g)
directory = '/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Python_Master_Thesis'
plt.savefig(directory + "/Electron_Hole_Gen.pdf", format='pdf', dpi=1200)
