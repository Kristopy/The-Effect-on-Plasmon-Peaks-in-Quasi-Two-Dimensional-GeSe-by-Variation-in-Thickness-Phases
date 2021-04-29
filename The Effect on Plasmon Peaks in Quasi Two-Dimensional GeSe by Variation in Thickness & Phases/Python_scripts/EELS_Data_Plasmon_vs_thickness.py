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
from matplotlib.pyplot import cm
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})
def spectrum_EELS(intensity):
    spectrum_pure = intensity[0:len(intensity)//3]
    background = intensity[len(intensity)//3: 2*len(intensity)//3]
    spectrum_no_background = intensity[2*len(intensity)//3: 3*len(intensity)//3]
    return spectrum_pure, background, spectrum_no_background

def plasmon_vs_thickness(Loc, file, Total_length,Total_length_thickness, num_gauss ,Gaussian, range_index_gaus, peak, sigma):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
    directory = folder + Loc
    directory1 = folder + 'GeSe_flake/Relative_thickness/EELS_REL_thickness_low_loss_GeSe_flake_ROI_2_POS_1.txt'
    data_folder = Path(directory)
    data_folder_thickness = Path(directory1)
    plt.rcParams["figure.figsize"] = (15,5)
    fig, ax1 = plt.subplots()

    #Inset axis
    left, bottom, width, height = [0.05, 0.5, 0.25, 0.4]
    ax2 = fig.add_axes([left, bottom, width, height])

    n = 14
    colors = np.flipud(plt.cm.copper(np.linspace(0,1,n)))#viridis
    #colors = ['darkred','lightcoral', 'sandybrown', 'goldenrod','lightskyblue', 'lightpink', 'black']

    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename.endswith(file):
            if i >= 0 and i <= 13 and i != 7 and i != 11:
                f = open(data_folder / filename, 'r').read()
                pure_list = f.split()
                tot_intensity = [float(i) for i in pure_list]
                spectrum_pure, background, spectrum_no_background = spectrum_EELS(tot_intensity)
                Energy_eV = np.linspace(0,Total_length, len(spectrum_pure)) #len of spectrum_EELS should be the same
                norm = list(minmax_scale(np.array(spectrum_no_background))) #Choose which data_set one want, background, EELS

                final_filename = filename.replace("_", "\_") #For adding to labels or legend: f'{final_filename[33:-4]}' for thickness
                print(final_filename)

                #Collecting average thickness values of 10x10nm and 1.5x10nm data files. ALSO included ROI_1_POS_1 for 1.3x17nm inset - stored in average_thickness_15x10nm when correct thickness file is used.

                average_thickness_10x10nm, average_thickness_15x10nm = pure_thickenss('GeSe_flake/Relative_thickness', 'EELS_REL_thickness_low_loss_GeSe_flake_ROI_2_POS_1.txt',Total_length_thickness)

                #average_thickness_10x10nm, average_thickness_15x10nm = pure_thickenss('GeSe_flake/Relative_thickness', 'EELS_REL_thickness_low_loss_GeSe_flake_ROI_1_POS_1.txt',Total_length_thickness)


                xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func = gaussian_fit(Energy_eV,spectrum_no_background, num_gauss, range_index_gaus, peak, sigma)
                print(y_peak, center)

                ax1.plot(center[y_peak.index(max(y_peak))], y_peak[y_peak.index(max(y_peak))], marker = "o", markersize=2, color='k')
                ax2.plot(Energy_eV,Gaussian_func[y_peak.index(max(y_peak))], alpha = 0.4, color = 'grey')

                #ax2.plot(Energy_eV[:100],norm_new[:100], alpha = 0.4, c = colors[i])

                #direct   =  list(minmax_scale(np.array((Energy_eV[:100] - 1.7)**(1/2))))
                #inderect =  list(minmax_scale(np.array((Energy_eV[:100] - 1.7)**(3/2))))
                #ax2.plot(Energy_eV[:100],direct, alpha = 0.4, c = 'r')
                #ax2.plot(Energy_eV[:100],inderect, alpha = 0.4, c = 'g')
                #ax2.plot(Energy_eV[:100],(Energy_eV[:100] - 1.7)**(1/2), alpha = 0.4, c = 'g')
                '''
                - INSET:
                you can manipulate the inset plot by either chosing gaussians, normalized function by estimated gassians, or real data.
                Typically, when chosing normalized function by estimated gaussians or real data, one is interested in the tail, and use tauc plot to
                either see direct or inderect gap. As for now, the estimated functionf from a number of gaussians give a smoother function, which
                is easier to look at.

                -PLOTTING:
                Be carefull, When plotting make sure you have chosen correct folder[], followed by correct length and thickness which corrensponds to the folder file.
                Furthermore, There are two plots to chose from, depending on which file in folder you have chosen, average_thickness_15x10nm and average_thickness_10x10nm.
                By ticking of one of these you then choose which t/lambda to use from function def pure_thickness. Make sure you have ticked off correct
                function as well. Two to chose from, EELS_REL_thickness_low_loss_GeSe_flake_ROI_2_POS_1 and EELS_REL_thickness_low_loss_GeSe_flake_ROI_1_POS_1, depending on which
                folder[] you have chosen.

                Furthermore, you need to adjust the parameter n, which is the number of colors from nmap. It is rather simple.
                EELS_Spectrum_high_loss_10x10nm_GeSe_flake_ROI_2_POS_1 has 10 files, meaning n = 11 is sufficient.
                - n = 11
                - if i >= 0:
                EELS_Spectrum_high_loss_1.32x17.4nm_GeSe_flake_ROI_1_POS_1 has 10 files, meaning n = 11 is sufficient.
                - n = 11
                - if i >= 0:
                EELS_Spectrum_high_loss_1.5x10nm_GeSe_flake_ROI_2_POS_1 has 20 files, meaning n = 21 is sufficient, NOTE: Loop is cut of at 13, meaning n = 14 is good enough. for accurate colors.
                - n = 14
                - i >= 0 and i <= 13 and i != 7:
                '''
                ax1.plot(Energy_eV[:-1200], spectrum_no_background[:-1200]  , label = 't/$\\lambda$: {:.2f}, $E_p$ {:.2f}eV'.format(average_thickness_15x10nm[i], center[y_peak.index(max(y_peak))]), c = colors[i])
                #ax1.plot(Energy_eV[:-1200], spectrum_no_background[:-1200]  , label = 't/$\\lambda$: {:.2f}, $E_p$ {:.2f}eV'.format(average_thickness_10x10nm[i], center[y_peak.index(max(y_peak))]), c = colors[i])


                ax2.axes.yaxis.set_visible(False)
                ax1.set_xlabel("Energy [eV]",fontsize=15)
                ax1.set_ylabel('Intensity [a.u]',fontsize=17)
                ax1.set_yticks([])
                ax1.set_xticks([5,10,15,20,25,30,35,40,45,50])
                ax2.set_xticks([5,10,15,20,25,30,35])
                ax2.set_xlim([0, 30])
                #ax2.set_xlim([min(Energy_eV[:100]), max(Energy_eV[:100])])

                plt.title('', fontsize=17)
                ax1.legend(loc='upper right', prop={"size": 17},ncol=2,  frameon=False)
                ax1.set_xlim([min(Energy_eV[:-1200]), max(Energy_eV[:-1200])])
                ax1.set_ylim([min(spectrum_no_background), max(spectrum_no_background) + max(spectrum_no_background)/2 ])
                plt.tight_layout()
                plt.savefig(folder + 'TEST_1.5_10'+ ".pdf", format='pdf', dpi=1200)


def pure_thickenss(Loc, file, Total_length_thickness):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
    directory = folder + Loc
    data_folder = Path(directory)

    f = open(data_folder/file, 'r').read()

    pure_list_thickenss =f.split()
    tot_intensity_thickness = [float(i) for i in pure_list_thickenss]
    thickness_nm = list(np.linspace(0,Total_length_thickness, len(tot_intensity_thickness)))
    tot_intensity_thickness = np.flipud(tot_intensity_thickness)
    final_filename = file.replace("_", "\_") #For adding to labels or legend: f'{final_filename[33:-4]}' for thickness
    print('thickness: ', thickness_nm)
    average_thickness_10x10nm = []
    average_thickness_15x10nm = []

    for i in range(11):
        average_thickness_10x10nm.append(np.mean(tot_intensity_thickness[23*i: 23*(i+1)]))
        print('thickeness 1010:', np.mean(tot_intensity_thickness[23*i: 23*(i+1)]))
    for i in range(21):
        if file == 'EELS_REL_thickness_low_loss_GeSe_flake_ROI_1_POS_1.txt':
            average_thickness_15x10nm.append(np.mean(tot_intensity_thickness[4*i: 4*(i+1)]))
        if file == 'EELS_REL_thickness_low_loss_GeSe_flake_ROI_2_POS_1.txt':
            average_thickness_15x10nm.append(np.mean(tot_intensity_thickness[3*i: 3*(i+1)]))

    return average_thickness_10x10nm, average_thickness_15x10nm

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
    print('Total number of peaks: ', len(peaks))
    print('-----------------------------------')
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

        print('-------------------------')
        print('Gaussian peak nr: ', i+1)
        print('Center: ', x[peaks[i]])
        print('mu: ',mu)
        print('sigma: ', sigma)
        print('rms: ', rms)
        print('x-span: ', x[peak_start[0]], x[peak_end[0]])
        print('Rel Intensity: ', y_1[peaks[i]])


    xGaus = x[peaks[peak]-range_1:peaks[peak]+range_1]
    norm_new = y_1[peaks[peak]-range_1:peaks[peak]+range_1]

    return xGaus, y_1, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func

#nm
GeSe_flake_spectrum_lenght = [103.504, 103.551] #nm 103.551 for GeSe_flake_ROI_1_POS_1.... 103.504 for GeSe_flake_ROI_2_POS_1
Total_length_thickness = [179.82, 51.03] #51.03 for GeSe_flake_ROI_1_POS_1 .... 179.82 for GeSe_flake_ROI_2_POS_1
Folder = ['GeSe_flake/EELS_Spectrum_high_loss_1.5x10nm_GeSe_flake_ROI_2_POS_1', 'GeSe_flake/EELS_Spectrum_high_loss_10x10nm_GeSe_flake_ROI_2_POS_1', 'GeSe_flake/EELS_Spectrum_high_loss_1.32x17.4nm_GeSe_flake_ROI_1_POS_1','GeSe_flake/Relative_thickness']
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

plasmon_vs_thickness(Folder[0], '.txt',GeSe_flake_spectrum_lenght[1],Total_length_thickness[0],  200,False,100, 0, sigma)
