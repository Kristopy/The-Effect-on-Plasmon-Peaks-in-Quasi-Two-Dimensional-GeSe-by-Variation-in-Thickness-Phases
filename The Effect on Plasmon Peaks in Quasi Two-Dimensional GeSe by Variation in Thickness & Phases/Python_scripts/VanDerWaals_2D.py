import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})

def ljp(r, epsilon, sigma):

    return 48 * epsilon * np.power(sigma, 12) / np.power(r, 13) \
    - 24 * epsilon * np.power(sigma, 6) / np.power(r, 7)

r = np.linspace(3.5, 8, 100)
plt.plot(r, ljp(r, 0.0103, 3.3))
plt.xlabel('distance')
plt.ylabel('energy')
plt.show()
