
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib import rc, rcParams
import os
import itertools
from pathlib import Path
import sys
from lmfit.models import ExponentialModel, GaussianModel
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")

def Reading_files(lines):
    x = []
    y = []
    for i in range(len(lines)-1):
        if (i % 2 == 0 and i != 0):
            splitted_lines = lines[i].split(',')
            c = list(map(float, splitted_lines))
            x.append(c[0])  # First element x-values
            y.append(c[1])  # second element u(x)
    return x,y

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt


folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/CL-Analysis/GeSe_4_Data/"
directory = folder + "Loc_4"
data_folder = Path(directory)

for filename in sorted(os.listdir(directory)):
    if filename.endswith('.txt'):
        print(filename)
        f = open(data_folder / filename, 'r').read()
        xdata,ydata = np.array(Reading_files(f.split("\n")))


# generate simulated data
#np.random.seed(123)  # comment out if you want different data each time
#xdata = np.linspace(0, 10, 100)

#ydata_perfect = gauss(xdata, 20, 5, 6, 1)
#ydata = np.random.normal(ydata_perfect, 1, 100)

H, A, x0, sigma = gauss_fit(xdata, ydata)

FWHM = 2.35482 * sigma

print('The offset of the gaussian baseline is', H)
print('The center of the gaussian fit is', x0)
print('The sigma of the gaussian fit is', sigma)
print('The maximum intensity of the gaussian fit is', H + A)
print('The Amplitude of the gaussian fit is', A)
print('The FWHM of the gaussian fit is', FWHM)

'''
#Try to use GaussionModel built in, then approximate Sigma.

folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/CL-Analysis/GeSe_4_Data/"
directory = folder + "Loc_1"
data_folder = Path(directory)

for filename in sorted(os.listdir(directory)):
    if filename.endswith('.txt'):
        print(filename)
        f = open(data_folder / filename, 'r').read()
        xdata,ydata = np.array(Reading_files(f.split("\n")))



#Parameters for loc_3 2.txt
'''
ix1 = index_of(xdata, 400)
ix2 = index_of(xdata, 590)
ix3 = index_of(xdata, 610)
ix4 = index_of(xdata, 710)
'''

exp_mod = ExponentialModel(prefix='exp_')
gauss1 = GaussianModel(prefix='g1_')
gauss2 = GaussianModel(prefix='g2_')
gauss3 = GaussianModel(prefix='g3_')
gauss4 = GaussianModel(prefix='g4_')



def index_of(arrval, value):
    """return index of array *at or below* value """
    if value < min(arrval):
        return 0
    return max(np.where(arrval <= value)[0])


ix1 = index_of(xdata, 290)
ix2 = index_of(xdata, 430)
ix3 = index_of(xdata, 530)
ix4 = index_of(xdata, 610)
ix5 = index_of(xdata, 710)
print(ix1, ix2, ix3, ix4 )

pars1 = exp_mod.guess(ydata[:ix1], x=xdata[:ix1])
pars2 = gauss1.guess(ydata[ix1:ix2], x=xdata[ix1:ix2])
pars3 = gauss2.guess(ydata[ix2:ix3], x=xdata[ix2:ix3])
pars4 = gauss3.guess(ydata[ix3:ix4], x=xdata[ix3:ix4])
pars5 = gauss4.guess(ydata[ix4:ix5], x=xdata[ix4:ix5])


pars =pars1+ pars2 + pars3 + pars4 + pars5
mod =gauss1 + gauss2 + gauss3 + gauss4 + exp_mod

out = mod.fit(ydata, pars, x=xdata)

print(out.fit_report(min_correl=0.5))

comps = out.eval_components(x=xdata)

#plt.plot(x, out.init_fit, 'k--', label='initial fit')


plt.plot(xdata, ydata, 'ko', label='data', markersize=3)
#plt.plot(xdata, ydata_perfect, '-k', label='data (without_noise)')
#plt.plot(xdata, gauss(xdata, *gauss_fit(xdata, ydata)), '--r', label='fit')

plt.plot(xdata, out.best_fit, 'r-', label='best fit')
plt.plot(xdata, comps['g1_'], 'darkred', )
plt.plot(xdata, comps['g2_'], 'lightskyblue')
plt.plot(xdata, comps['g3_'], 'orange')
plt.plot(xdata, comps['g4_'], 'black')
#plt.plot(xdata, comps['exp_'])

plt.fill_between(xdata, comps['g1_'],color='darkred', alpha=0.2, label='Gaussian;1')
plt.fill_between(xdata, comps['g2_'],color='lightskyblue', alpha=0.2, label='Gaussian;2')
plt.fill_between(xdata, comps['g3_'],color='orange' , alpha=0.2, label='Gaussian;3')
plt.fill_between(xdata, comps['g4_'],color='black' , alpha=0.2, label='Gaussian;4')


plt.legend(loc='best')

directory1 = "/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/CL-Analysis/GeSe_4_Data"
plt.savefig(directory1 + "/TEST3.pdf", format='pdf', dpi=1200)
