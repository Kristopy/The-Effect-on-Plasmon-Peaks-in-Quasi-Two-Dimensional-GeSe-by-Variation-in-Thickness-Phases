import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
import os
import itertools
from pathlib import Path
import sys
from matplotlib import pyplot as plt
import math
from scipy.signal import find_peaks , peak_widths
from lmfit.models import ExponentialModel, GaussianModel
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")
from astroML.datasets import fetch_vega_spectrum
from astroML.sum_of_norms import sum_of_norms, norm
from sklearn.preprocessing import minmax_scale


rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})


def spectrum_EELS(intensity):
    spectrum_pure = intensity[0:len(intensity)//3]
    spectrum_no_background = intensity[len(intensity)//3: 2*len(intensity)//3]
    background = intensity[2*len(intensity)//3: 3*len(intensity)//3]
    return spectrum_pure, spectrum_no_background, background

def Position(file, Total_length, num_gauss ,Gaussian, model, range_index_gaus, peak):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
    directory = folder
    data_folder = Path(directory)
    j = 0
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            f = open(data_folder / filename, 'r').read()
            pure_list = f.split()
            tot_intensity = [float(i) for i in pure_list]
            spectrum_pure, spectrum_no_background, background = spectrum_EELS(tot_intensity)

            Energy_eV = np.linspace(0,10, len(tot_intensity[0:len(tot_intensity)//3]))

            #l_b_p = Lenght_peaks/(Number_of_peaks-1)
            colors = ['darkred', 'goldenrod','lightskyblue', 'lightpink', 'black']
            plt.rcParams["figure.figsize"] = (15,5)
            if Gaussian == False and model in filename and model == 'Spectrum':
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    print(filename[:-4])
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    fig, ax1 = plt.subplots()

                    #Inset axis
                    #left, bottom, width, height = [0.3, 0.6, 0.3, 0.3]
                    #ax2 = fig.add_axes([left, bottom, width, height])

                    Total_length = [103.551, 103.03, 103.874, 103.504, 103.504] #nm
                    Lenght_x = np.linspace(0,Total_length[j], len(spectrum_pure))
                    norm = list(minmax_scale(np.array(background)))


                    lines = gaussian_fit(Lenght_x,norm, num_gauss, range_index_gaus, ax1, peak)
                    ax1.plot(Lenght_x,norm, linestyle='-',linewidth=0.7)

                    plt.savefig(folder + filename[:-4] + '_Gaussian' + ".pdf", format='pdf', dpi=1200)


def gaussian_fit(x,y, num_gauss, range_1, ax, peak):

    x =np.array(x)
    y =np.array(y)

    for n_gaussians in range(1,num_gauss):
        # compute the best-fit linear combination
            w_best, rms, locs, widths = sum_of_norms(x, y, n_gaussians,
                                                     spacing='linear',
                                                    full_output=True)
            norms = w_best * norm(x[:, None], locs, widths)

            if rms <= 0.00000001:
                print('Number of Gaussians: ', n_gaussians)
                break

    y_1 = norms.sum(1)
    grad = np.gradient(y_1)
    ax.plot(x, y_1)

    peaks, _ = find_peaks(y_1[:-800],height=0.01, width=1) #Set 0.19 for GeSe_flake
    results_half = peak_widths(y_1, peaks, rel_height=0.005)

    print('Total number of peaks: ', len(peaks))
    print('-----------------------------------')
    lines_gaus = []
    for i in range(len(peaks)):

        peak_start = np.where(x >= abs(x[peaks[i]] - float(results_half[0][i])/2))[0]
        peak_end   = np.where(x >= abs(x[peaks[i]] + float(results_half[0][i])/2))[0]

        ax.plot(x[peaks[i]], y[peaks[i]], 'o', markersize=3, color='black')

        mu = x[peaks[i]]
        sigma = np.std(x[peak_start[0]: peak_end[0]], ddof=1)
        Gaus = 0.5*y_1[peaks[i]] * np.exp(-(x - mu)**2 / (2 * sigma**2))
        gaus_shift = [i-min(Gaus) for i in Gaus]

        lines = ax.fill_between(x, gaus_shift , alpha=0.4)
        lines_gaus.append(lines)
    return lines_gaus


Total_length_GeSe_Needle_spectrum = [103.479, 103.628, 103.704, 103.828] #nm

Position('.txt',Total_length_GeSe_Needle_spectrum,  200,False, 'Spectrum',10, 0 )
