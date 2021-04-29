import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})


t = np.linspace(0,1,10e4)

def Recombination(t):
    delta_n = 1

    n = delta_n*np.exp(-(t)/0.2)

    return n


n = Recombination(t)

plt.figure()

plt.plot(t, n, linestyle='-',linewidth=0.7, color='darkred')

plt.legend(loc='upper right', prop={"size": 13}, frameon=False)
plt.xlabel("Recombination time [ns]",fontsize=15)
plt.ylabel('$\\delta n(t)$, $\\delta p(t)$',fontsize=15)
plt.autoscale(enable=True, axis='x', tight=True)
#plt.xticks([])
plt.yticks([])
#plt.title("Height profile of GeSe; Position I", fontsize=17)


directory = '/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Python_Master_Thesis'
plt.savefig(directory + "/Recombination.eps", format='eps', dpi=1200)
