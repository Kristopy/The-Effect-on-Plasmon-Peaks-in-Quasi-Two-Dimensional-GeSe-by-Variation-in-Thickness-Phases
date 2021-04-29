import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
import os
import itertools
from pathlib import Path
import sys
from matplotlib import pyplot as plt
import math
from scipy.signal import find_peaks, peak_widths
from lmfit.models import ExponentialModel, GaussianModel
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")
from astroML.datasets import fetch_vega_spectrum
from astroML.sum_of_norms import sum_of_norms, norm
from sklearn.preprocessing import minmax_scale
from sklearn.linear_model import LinearRegression
from Atomic_percents_EDX import Atomic_percent

from matplotlib.pyplot import cm

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})
def spectrum_EELS(intensity):
    spectrum_pure = intensity[0:len(intensity)//3]
    background = intensity[len(intensity)//3: 2*len(intensity)//3]
    spectrum_no_background = intensity[2*len(intensity)//3: 3*len(intensity)//3]
    return spectrum_pure, background, spectrum_no_background

def gaussian_fit(x,y, num_gauss, range_1, peak, sigma):

    x =np.array(x)

    y = np.array(y)
    y1 =minmax_scale(np.array(y))

    for n_gaussians in range(1,num_gauss):
        # compute the best-fit linear combination
            w_best, rms, locs, widths = sum_of_norms(x, y1, n_gaussians,
                                                     spacing='linear',
                                                    full_output=True)
            norms = w_best * norm(x[:, None], locs, widths)

            if rms <= 10e-10:
                print('Number of Gaussians: ', n_gaussians)
                break

    y_1 = norms.sum(1)

    # compute second derivative
    #second_derivative = np.gradient(np.gradient(y_1[:-1100]))

    # find switching points
    #infls = np.where(np.diff(np.sign(second_derivative)))[0]
    #print(np.sign(second_derivative))
    #ax.plot([x[infls], x[infls]], [0,1])

    peaks, _ = find_peaks(y_1[:-930],height=0.01, width=1)
    results_half = peak_widths(y_1, peaks, rel_height=0.001)
    #peaks = np.delete(peaks, 2) #Deleting peaks by index
    #print('Total number of peaks: ', len(peaks))
    #print('-----------------------------------')
    #peaks = peaks[:3] + peaks[3+1 :]

    center = []
    mu_1 = []
    sigma_1 = []
    y_peak = []
    Gaussian_func = []
    for i in range(len(peaks)):
        peak_start = np.where(x >= abs(x[peaks[i]] - float(results_half[0][i])/2))[0]
        peak_end   = np.where(x >= abs(x[peaks[i]] + float(results_half[0][i])/2))[0]

        mu = x[peaks[i]]

        #sigma = np.std(x[peak_start[0]: peak_end[0]], ddof=1)
        Gaus = 0.5*y[peaks[i]] * np.exp(-(x - mu)**2 / (2 * sigma[i]**2))


        Gaussian_func.append(Gaus)
        center.append(x[peaks[i]])
        mu_1.append(mu)
        sigma_1.append(sigma[i])
        y_peak.append(y[peaks[i]])

        #print('-------------------------')
        #print('Gaussian peak nr: ', i+1)
        #print('Center: ', x[peaks[i]])
        #print('mu: ',mu)
        #print('sigma: ', sigma)
        #print('rms: ', rms)
        #print('x-span: ', x[peak_start[0]], x[peak_end[0]])
        #print('Rel Intensity: ', y_1[peaks[i]])


    xGaus = x[peaks[peak]-range_1:peaks[peak]+range_1]
    norm_new = y_1[peaks[peak]-range_1:peaks[peak]+range_1]

    return xGaus, y_1, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func


def plasmon_vs_thickness(Loc, file, Total_length, num_gauss ,Gaussian, range_index_gaus, peak, sigma):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
    directory = folder + Loc
    data_folder = Path(directory)
    j = 0
    lab1 = []
    print(directory)
    plt.rcParams["figure.figsize"] = (15,5)


    #Inset axis
    #left, bottom, width, height = [0.05, 0.5, 0.25, 0.4]
    #plt = fig.add_axes([left, bottom, width, height])
    n = 6
    colors = np.flipud(plt.cm.copper(np.linspace(0,1,n)))#viridis

    #colors = ['darkred','lightcoral', 'sandybrown', 'goldenrod','lightskyblue', 'lightpink', 'black']
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename.endswith(file):
            f = open(data_folder / filename, 'r').read()
            pure_list = f.split()
            tot_intensity = [float(i) for i in pure_list]
            spectrum_pure, background, spectrum_no_background = spectrum_EELS(tot_intensity)
            Energy_eV = np.linspace(0,Total_length[j], len(spectrum_pure)) #len of spectrum_EELS should be the same
            norm = list(minmax_scale(np.array(spectrum_no_background))) #Choose which data_set one want, background, EELS
            final_filename = filename.replace("_", "\_") #For adding to labels or legend: f'{final_filename[33:-4]}' for thickness

            final_filename = final_filename[33:-4]
            print(final_filename)
            xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func = gaussian_fit(Energy_eV,spectrum_no_background, num_gauss, range_index_gaus, peak, sigma)

            #Linear regression:

            Energy_eV_1 = Energy_eV[40:55]
            norm_1 = norm_new[40 :55]
            lin_reg_Energy_eV = np.array(Energy_eV_1).reshape((-1, 1))

            model = LinearRegression().fit(lin_reg_Energy_eV, norm_1)

            r_sq = model.score(lin_reg_Energy_eV, norm_1)
            print('coefficient of determination:', r_sq)
            print('intercept:', model.intercept_)
            print('slope:', model.coef_)
            Eg = -(model.intercept_)/model.coef_
            print('band gap: ', Eg)

            dict_atoms = Atomic_percent('.txt', '/NEW_TEST.pdf')

            if final_filename in dict_atoms.keys():
                print('dict_atoms: ', dict_atoms[final_filename])

                lab1.append('$E_p$: {:.2f} $E_g$:  {:.2f}  $Ge_{{{:.2f}}}$$Se_{{{:.2f}}}$'.format(center[y_peak.index(max(y_peak))], Eg[0], dict_atoms[final_filename][0]/100, dict_atoms[final_filename][1]/100))
            else:
                print('OK')
                lab1.append('{} : $E_p$: {:.2f} $E_g$:  {:.2f}'.format(final_filename, center[y_peak.index(max(y_peak))], Eg[0]))

            plt.plot(center[y_peak.index(max(y_peak))], norm[y_peak.index(max(y_peak))], marker = "o", markersize=2, color='k')
            plt.plot(Energy_eV[:100], norm_new[:100], c = colors[j], label=lab1[j])
            plt.plot(Energy_eV_1, model.intercept_+ model.coef_*Energy_eV_1, 'r', linewidth=0.7, linestyle='--')
            plt.legend(loc='upper left', frameon=False)


            #plt.plot(Energy_eV[:100], norm[:100]  , label = ' - center {:.2f}'.format(center[y_peak.index(max(y_peak))]), c = colors[0])
            plt.xlabel("Energy [eV]",fontsize=15)
            plt.ylabel('Intensity [a.u]',fontsize=17)
            plt.title('', fontsize=17)
            plt.legend(loc='upper left', prop={"size": 17},ncol=2,  frameon=False)
            plt.yticks([])
            plt.xlim([min(Energy_eV[:80]), max(Energy_eV[:80])])
            plt.ylim([min(norm[:60]), max(norm[:60]) + max(norm[:60])/2 ])
            plt.tight_layout()
            plt.savefig(folder + 'TEST'+ ".pdf", format='pdf', dpi=1200)
            j += 1

#nm
GeSe_flake_spectrum_lenght = [103.504, 103.551] #nm 103.551 for GeSe_flake_ROI_1_POS_1.... 103.504 for GeSe_flake_ROI_2_POS_1

Total_length_GeSe_flake_spectrum = [103.551, 103.03, 103.874, 103.504, 103.504] #nm
Total_length_GeSe_acu_spectrum = [103.905, 103.826, 103.807, 103.906] #nm
Total_length_GeSe_Needle_spectrum = [103.479, 103.628, 103.704, 103.828] #nm

Folder = ['GeSe_flake/EELS_Spectrum_high_loss_1.5x10nm_GeSe_flake_ROI_2_POS_1', 'GeSe_flake/EELS_Spectrum_high_loss_10x10nm_GeSe_flake_ROI_2_POS_1', 'GeSe_flake/EELS_Spectrum_high_loss_1.32x17.4nm_GeSe_flake_ROI_1_POS_1']
folders = ['GeSe_flake','GeSe_accumulation','GeSe_Needle']

sigma = [4,2,2,2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]

#Defining arguments in position for later use:
# 1: Folder to execute script
# 2: ending of file for execution: for ex - if file ends with .txt, placing '.txt' in arguments will make the code run thorugh all .txt file in that given folder.
# 3: Set total length for what you want to plot.
# 4: Number of Gaussians for estimating the function from .txt for better estimate of peaks.
# 5: True or False argument for: True means that we dont esitmate gaussians and peaks. False mean we do execute Gaussian model and find peaks.
# 6: Two posibilities so far: 'Spectrum' finds all .txt which has 'Spectrum' as a word - this contains only High loss spectrums. Or 'thickness' which plots thickness profile along a direction
# Alse added to nr.6 is 'peak_table_plot' - when nr.5 is set on False and nr.6 ia 'peak_table_plot' we plot data of all plasmon peaks and compare.
# 7: Setting number of indexes from peak position to be displayed in inset.
# 8: inset peak- insert an index - 0 ... peak.


#plasmon_vs_thickness(Folder[0], '.txt',GeSe_flake_spectrum_lenght[0],200,False,100, 0, sigma)
plasmon_vs_thickness(folders[0], '.txt',Total_length_GeSe_flake_spectrum,200,False,100, 0, sigma)
