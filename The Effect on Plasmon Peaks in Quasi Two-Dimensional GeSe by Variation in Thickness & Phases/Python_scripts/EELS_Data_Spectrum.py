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
from Atomic_percents_EDX import Atomic_percent

#rc('text', usetex=True)
#rc('font', **{'family': 'serif', 'serif': ['Random']})


def spectrum_EELS(intensity):
    spectrum_pure = intensity[0:len(intensity)//3]
    background = intensity[len(intensity)//3: 2*len(intensity)//3]
    spectrum_no_background = intensity[2*len(intensity)//3: 3*len(intensity)//3]
    return spectrum_pure, background, spectrum_no_background

def EELS_Data_Spectrum(Loc, file, Total_length, num_gauss ,Gaussian, model, range_index_gaus, peak, sigma):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
    directory = folder + Loc
    data_folder = Path(directory)
    j = 0
    lines_comb1 = []
    lines_comb2 = []
    lines_comb3 = []
    lab1 = []
    lab2 = []
    lab3 = []
    plt.rcParams["figure.figsize"] = (15,5)
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            f = open(data_folder / filename, 'r').read()
            pure_list = f.split()
            tot_intensity = [float(i) for i in pure_list]
            spectrum_pure, background, spectrum_no_background = spectrum_EELS(tot_intensity)
            Energy_eV = np.linspace(0,Total_length[j], len(spectrum_pure)) #len of spectrum_EELS should be the same
            Energy_eV_thickness = np.linspace(0,Total_length[j], len(tot_intensity)) #x-data for thickness measurments
            norm = list(minmax_scale(np.array(spectrum_no_background))) #Choose which data_set one want, background, EELS
            colors = ['darkred', 'goldenrod','lightskyblue', 'lightpink', 'black']


            final_filename = filename.replace("_", "\_") #For adding to labels or legend: f'{final_filename[33:-4]}' for thickness

            if Gaussian == True:
                if model in filename and model == 'thickness':
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    print(filename[:-4])
                    print('-------------------------------------------')
                    print('-------------------------------------------')

                    plt.rcParams["figure.figsize"] = (15,5)
                    plt.figure()

                    plt.fill_between(Energy_eV_thickness, tot_intensity  , alpha=0.4, color = colors[j], label = f'{final_filename[33:-4]}')#color = colors[j-2]

                    print('thickenss')
                    plt.xlabel("Spectral image length [nm]",fontsize=15)
                    plt.ylabel('$\\frac{t}{\\lambda}$',fontsize=17)
                    #plt.title('', fontsize=17)
                    plt.legend(loc='upper right', prop={"size": 17}, frameon=False)
                    plt.xlim([min(Energy_eV_thickness), max(Energy_eV_thickness)])
                    plt.ylim([min(tot_intensity), max(tot_intensity) + max(tot_intensity)/10 ])
                    plt.tight_layout()
                    #plt.savefig(folder + Loc + '/Figures/' + filename[:-4] + ".pdf", format='pdf', dpi=1200)
                    plt.savefig(folder + Loc + '/Figures/' + 'TEST' + ".pdf", format='pdf', dpi=1200)
                    j += 1

                if model in filename and model == 'Spectrum':
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    print(filename[:-4], Loc)
                    print('-------------------------------------------')
                    print('-------------------------------------------')


                    dict_atoms = Atomic_percent('.txt', '/NEW_TEST.pdf')

                    xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func = gaussian_fit(Energy_eV,norm, num_gauss, range_index_gaus, peak, sigma)

                    if Loc == 'GeSe_flake':

                        if np.where(Total_length.index(Total_length[j]) == 0)[0] == 0:
                            plt.figure()
                        lab1.append('{}: {:.2f}eV'.format(final_filename[33:-4], center[y_peak.index(max(y_peak))]))
                        lines_comb1 += plt.plot(Energy_eV, spectrum_no_background  ,color = colors[j])
                        leg1 = plt.legend(lines_comb1,lab1, loc='upper right', frameon=False)
                    if Loc == 'GeSe_accumulation':
                        if np.where(Total_length.index(Total_length[j]) == 0)[0] == 0:
                            plt.figure()

                        lab2.append('{}: {:.2f}eV'.format(final_filename[33:-4], center[y_peak.index(max(y_peak))]))
                        lines_comb2 += plt.plot(Energy_eV, spectrum_no_background  , color = colors[j] )
                        leg2 = plt.legend(lines_comb2,lab2, loc='upper right', frameon=False)
                    if Loc == 'GeSe_Needle':
                        if np.where(Total_length.index(Total_length[j]) == 0)[0] == 0:
                            plt.figure()
                        lab3.append('{}: {:.2f}eV'.format(final_filename[33:-4], center[y_peak.index(max(y_peak))]))
                        lines_comb3 += plt.plot(Energy_eV, spectrum_no_background  ,  color = colors[j] )
                        leg3 = plt.legend(lines_comb3,lab3, loc='upper right', frameon=False)

                    #n = 9
                    #colors = np.flipud(plt.cm.viridis(np.linspace(0,1,n)))#viridis
                    #lab1.append('{}: {:.2f}eV'.format(final_filename[33:-4], center[y_peak.index(max(y_peak))]))
                    #plt.plot(Energy_eV, spectrum_no_background  ,color = colors[j])
                    plt.legend(lab1, loc='upper right', frameon=False)

                    plt.xlabel("Energy [eV]",fontsize=15)
                    plt.ylabel('Intensity [a.u]',fontsize=15)
                    plt.title('', fontsize=17)
                    plt.xlim([min(Energy_eV), max(Energy_eV)])
                    #plt.ylim([min(norm), max(norm) + max(norm)/10 ])
                    plt.tight_layout()
                    plt.savefig(folder + 'Combined' + ".pdf", format='pdf', dpi=1200)
                    j += 1

            if Gaussian == False and model in filename and model == 'Spectrum':
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    print(filename[:-4])
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    plt.rcParams["figure.figsize"] = (15,5)

                    fig, ax1 = plt.subplots()

                    #Inset axis
                    left, bottom, width, height = [0.5, 0.6, 0.3, 0.3]
                    ax2 = fig.add_axes([left, bottom, width, height])

                    xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func = gaussian_fit(Energy_eV,norm, num_gauss, range_index_gaus, peak, sigma)
                    label_gaus = []

                    lines_gaus = []
                    for i in range(len(center)):
                        lines = ax1.fill_between(Energy_eV, Gaussian_func[i] , alpha=0.4, color = colors[i])
                        lines_gaus.append(lines)
                        label_gaus.append('\\bf{Gaussian peak} %.0f:\n Center $\\mu$: %.2f [eV] \n$\\sigma : $ %.2f' % (i+1 ,center[i], sigma_1[i]))
                        #label_gaus.append('Center : {:.2f}eV\n  $\\mu    :$ {:>16}\n  $\\sigma : $ {:>16}'.format(center[i],mu_2[i],sigma_2[i]))


                    #By rewriting bbox_to_anchor you may move label position.
                    leg2 = ax1.legend(lines_gaus, label_gaus,
                              loc='upper right', frameon=False, fontsize=14, bbox_to_anchor=(0.975 ,0.95))
                    print(leg2)
                    #Basic plot of length
                    ax1.plot(Energy_eV,norm, linestyle='-',linewidth=0.7, color = colors[j])
                    #Plot inset for xGaus which is the width of from centrum of peak, distance controlled by range_index_gaus and norm_new which is the intensity values after Gaussian estimation.
                    ax2.plot(xGaus,norm_new, linestyle='-',linewidth=0.7, color = colors[j])

                    #Setting lines from inset box to peak

                    if peak == 0:
                        ax1.plot([48,center[peak]], [0.82  , y_peak[peak]+y_peak[peak]/20], linestyle='--', linewidth=0.3, color='black')
                        ax1.errorbar( [min(xGaus),max(xGaus)]   , [y_peak[peak]+y_peak[peak]/20, y_peak[peak]+y_peak[peak]/20], yerr=0.02, linewidth=0.3, color='black')
                    else:
                        ax1.plot([65,center[peak]], [0.65, y_peak[peak]+y_peak[peak]/20], linestyle='--', linewidth=0.3, color='black')
                        ax1.errorbar( [min(xGaus),max(xGaus)]   , [y_peak[peak]+y_peak[peak]/20, y_peak[peak]+y_peak[peak]/20], yerr=0.02, linewidth=0.3, color='black')

                    ax2.axes.yaxis.set_visible(False)
                    ax2.axes.xaxis.set_visible(False)
                    ax1.set_title('Number of Gaussians: %.0f -- rms error: %.2E' % (n_gaussians, rms), fontsize= 17)
                    ax1.set_xlabel("Energy [eV]",fontsize=15)
                    ax1.set_ylabel('Norm. Intensity [a.u]',fontsize=15)
                    #ax1.title('', fontsize=17)
                    ax1.set_xlim([min(Energy_eV), max(Energy_eV)])
                    ax1.set_ylim([min(norm), max(norm) + max(norm)/10 ])

                    ax2.set_xlim([min(xGaus), max(xGaus)])
                    ax2.set_ylim([min(norm_new), max(norm_new) + max(norm_new)/100  ])
                    plt.tight_layout()
                    #plt.savefig(folder + Loc + '/Figures/' + filename[:-4] + '_Gaussian' + ".pdf", format='pdf', dpi=1200, bbox_extra_artists=(leg2,), bbox_inches='tight')
                    #plt.savefig(folder + Loc + '/Figures/' + filename[:-4] + '_Gaussian' + ".pdf", format='pdf', dpi=1200, bbox_extra_artists=(leg2,), bbox_inches='tight')
                    j += 1

            if Gaussian == False and 'thickness' not in filename and model == 'peak_table_plot':
                    #print('-------------------------------------------')
                    #print('-------------------------------------------')
                    #print(filename[:-4])
                    #print('-------------------------------------------')
                    #print('-------------------------------------------')
                    plt.rcParams["figure.figsize"] = (15,5)
                    fig, ax1 = plt.subplots()
                    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
                    ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

                    xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func, y_1 = gaussian_fit(Energy_eV,norm, num_gauss, range_index_gaus, peak, sigma)
                    lines_peaks_plot.append(center[sigma_1.index(max(sigma_1))])
                    list_filenames.append(final_filename)

                    count_num_files = [i+1 for i in range(len(lines_peaks_plot))]


                    lab1 = []
                    lab2 = []
                    lab3 = []
                    lines_1 = []
                    lines_2 = []
                    lines_3 = []
                    ax1.plot(count_num_files[0:5],lines_peaks_plot[0:5], linestyle='--',color='k', linewidth=0.7)
                    ax1.plot(count_num_files[9:13],lines_peaks_plot[9:13], linestyle='--',color='k', linewidth=0.7)
                    ax1.plot(count_num_files[5:9],lines_peaks_plot[5:9], linestyle='--',color='k', linewidth=0.7)

                    dict_atoms = Atomic_percent('.txt', '/NEW_TEST.pdf')


                    #pl.text(count_num_files[j]- 0.1,lines_peaks_plot[j]-lines_peaks_plot[j]/30,'$Ge_{%.2f}Se_{%.2f}$' % (dict_atoms[final_filename[33:-4]][0]/100, dict_atoms[final_filename[33:-4]][1]/100))


                    if len(count_num_files) == 13:
                        for i in range(len(count_num_files)):
                            print(list_filenames[i])

                            if list_filenames[i][33:-4] in dict_atoms.keys():
                                print(dict_atoms[list_filenames[i][33:-4]])

                                ax1.annotate('$Ge_{%.2f}Se_{%.2f}$' % (dict_atoms[list_filenames[i][33:-4]][0]/100, dict_atoms[list_filenames[i][33:-4]][1]/100),
                                            xy=(count_num_files[i], lines_peaks_plot[i]),
                                            xytext=(0, 10),  # 3 points vertical offset
                                            textcoords="offset points",
                                            ha='center', va='bottom', size=14)

                            markers=[">", "^","s", "h", "+", "o", "v"]

                            if i <= 4:
                                if i == 4:
                                    lab1.append('{} : {:{width}.{prec}f}eV'.format(list_filenames[i][31:-4], lines_peaks_plot[i],width = 10,prec=2))
                                else:
                                    lab1.append('{} : {:{width}.{prec}f}eV'.format(list_filenames[i][33:-4], lines_peaks_plot[i],width = 10,prec=2))
                                lines_1 += ax1.plot(count_num_files[i],lines_peaks_plot[i],linestyle='--',color='k', linewidth=0.7 ,marker=markers[1],  mfc=colors[i],mec=colors[i])
                                leg1 =ax1.legend(lines_1,lab1, loc='upper left', frameon=False, prop={'size': 12})

                            elif 5 <= i <= 8:
                                if i == 8:
                                    lab2.append('{} : {:{width}.{prec}f}eV'.format(list_filenames[i][31:-4], lines_peaks_plot[i],width = 10, prec=2))
                                else:
                                    lab2.append('{} : {:{width}.{prec}f}eV'.format(list_filenames[i][33:-4], lines_peaks_plot[i],width = 10,prec=2))
                                lines_2 += ax1.plot(count_num_files[i],lines_peaks_plot[i],linestyle='--',color='k', linewidth=0.7 ,marker=markers[2] , mfc=colors[i-5], mec=colors[i-5])
                                leg2 =ax2.legend(lines_2,lab2, loc='upper center', frameon=False, prop={'size': 12})
                            else:
                                lab3.append('{} : {:{width}.{prec}f}eV'.format(list_filenames[i][33:-4], lines_peaks_plot[i],width = 10,prec=2))
                                lines_3 += ax1.plot(count_num_files[i],lines_peaks_plot[i],linestyle='--',color='k', linewidth=0.7,marker=markers[5],mfc=colors[i-9], mec=colors[i-9])
                                leg3 =ax3.legend(lines_3,lab3, loc='upper right', frameon=False, bbox_to_anchor=(0.925 ,0.25), prop={'size': 12})

                    ax2.set_axis_off()
                    ax3.set_axis_off()

                    ax1.plot([5.5,5.5],[min(lines_peaks_plot)-1, max(lines_peaks_plot)+1], linestyle='--', linewidth='0.7', color='k')
                    ax1.plot([9.5,9.5],[min(lines_peaks_plot)-1, max(lines_peaks_plot)+1], linestyle='--', linewidth='0.7', color='k')
                    ax1.set_ylabel('Plasmon peak energy [eV]',fontsize=17)
                    #ax1.set_title('Plasmon peaks for different ROI as a function of ', fontsize=17)
                    ax1.set_ylim([min(lines_peaks_plot)-1, max(lines_peaks_plot)+1])
                    ax1.set_xticks([])
                    #plt.ylim([min(norm), max(norm) + max(norm)/10 ])
                    plt.tight_layout()
                    plt.savefig(folder + 'HEI' + ".pdf", format='pdf', dpi=1200, bbox_inches='tight')

                    j += 1

def gaussian_fit(x,y, num_gauss, range_1, peak, sigma):

    x =np.array(x)
    y =np.array(y)

    for n_gaussians in range(1,num_gauss):
        # compute the best-fit linear combination
            w_best, rms, locs, widths = sum_of_norms(x, y, n_gaussians,
                                                     spacing='linear',
                                                    full_output=True)
            norms = w_best * norm(x[:, None], locs, widths)

            if rms <= 10e-10:
                print('Number of Gaussians: ', n_gaussians)
                break

    y_1 = norms.sum(1)

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

        sigma = np.std(x[peak_start[0]: peak_end[0]], ddof=1)
        Gaus = 0.5*y_1[peaks[i]] * np.exp(-(x - mu)**2 / (2 * sigma**2))


        Gaussian_func.append(Gaus)
        center.append(x[peaks[i]])
        mu_1.append(mu)
        sigma_1.append(sigma)
        y_peak.append(y_1[peaks[i]])
        '''
        print('-------------------------')
        print('Gaussian peak nr: ', i+1)
        print('Center: ', x[peaks[i]])
        print('mu: ',mu)
        print('sigma: ', sigma)
        print('x-span: ', x[peak_start[0]], x[peak_end[0]])
        print('Rel Intensity: ', y_1[peaks[i]])
        '''

    xGaus = x[peaks[peak]-range_1:peaks[peak]+range_1]
    norm_new = y_1[peaks[peak]-range_1:peaks[peak]+range_1]

    return xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func, y_1

'''
New peak table for single plot, insider EELS_Data_Spectrum you can plot in for loop, meaning you will
get all three sites plotted out. For def peak_table_plot you can choose one site at a time, meaning you could Choose
GeSe_flake, and only plot this.  -
'''
def peak_table_plots(Loc, file, Total_length, num_gauss ,Gaussian, model, range_index_gaus, peak, sigma):
        folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
        directory = folder + Loc
        data_folder = Path(directory)
        j = 0
        center_peak = []
        lines_1 = []
        text_1 = []
        lab1 = []
        leg1 = []
        lines_stp = []
        plt.rcParams["figure.figsize"] = (15,5)
        fig, ax1 = plt.subplots()

        for i, filename in enumerate(sorted(os.listdir(directory))):
            if filename.endswith(file):
                f = open(data_folder / filename, 'r').read()
                pure_list = f.split()
                tot_intensity = [float(i) for i in pure_list]
                spectrum_pure, background, spectrum_no_background = spectrum_EELS(tot_intensity)
                Energy_eV = np.linspace(0,Total_length[j], len(spectrum_pure)) #len of spectrum_EELS should be the same
                Energy_eV_thickness = np.linspace(0,Total_length[j], len(tot_intensity)) #x-data for thickness measurments
                norm = list(minmax_scale(np.array(spectrum_no_background))) #Choose which data_set one want, background, EELS
                colors = ['darkred', 'goldenrod','lightskyblue', 'lightpink', 'black']
                markers=[">", "^","s", "h", "+", "o", "v"]
                final_filename = filename.replace("_", "\_") #For adding to labels or legend: f'{final_filename[33:-4]}' for thickness

                if 'thickness' not in filename and model == 'peak_table_plot':

                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    print(filename[28:-4])
                    print('-------------------------------------------')
                    print('-------------------------------------------')
                    xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func = gaussian_fit(Energy_eV,norm, num_gauss, range_index_gaus, peak, sigma)
                    center_peak.append(center[sigma_1.index(max(sigma_1))])
                    lines_stp.append(j)
                    dict_atoms = Atomic_percent('.txt', '/NEW_TEST.pdf')
                    if final_filename[33:-4] in dict_atoms.keys():
                        print(dict_atoms[final_filename[33:-4]])

                        ax1.annotate('$Ge_{%.2f}Se_{%.2f}$' % (dict_atoms[final_filename[33:-4]][0]/100, dict_atoms[final_filename[33:-4]][1]/100),
                                    xy=(lines_stp[j], center_peak[j]),
                                    xytext=(0, 10),  # 3 points vertical offset
                                    textcoords="offset points",
                                    ha='center', va='bottom')

                        #pl.text(lines_stp[j]- 0.1,center_peak[j]-center_peak[j]/30,'$Ge_{%.2f}Se_{%.2f}$' % (dict_atoms[final_filename[33:-4]][0]/100, dict_atoms[final_filename[33:-4]][1]/100))


                    lab1.append('{} : {:{width}.{prec}f}eV'.format(final_filename[33:-4], center_peak[j],width = 10,prec=2))

                    lines_1 += ax1.plot(lines_stp[j],center_peak[j],linestyle='--',color='k', linewidth=0.7, marker=markers[j],markersize=8,  mfc=colors[j],mec=colors[j])
                    ax1.plot(lines_stp,center_peak,linestyle='--',color='k', linewidth=0.7)
                    leg1 =ax1.legend(lines_1,lab1, loc='upper left', frameon=False)
                    final_loc =  Loc.replace("_", "\_")
                    ax1.set_ylabel('Plasmon peak energy [eV]',fontsize=15)
                    ax1.set_xticks([])
                    ax1.set_title('Plasmon peak energies for %s' % (final_loc), fontsize=17)
                    ax1.set_ylim([min(center_peak)-1, max(center_peak)+1])
                    #ax1.set_ylim([min(center_peak)-1, max(lines_peaks_plot)+1])
                    plt.tight_layout()
                    plt.savefig(folder + '_Peaks_1' + ".pdf", format='pdf', dpi=1200, bbox_inches='tight')
                    j += 1

def total_EELS_data_spectrum(Loc, file, Total_length, num_gauss, range_index_gaus, peak, sigma):
        count = 0
        plt.rcParams["figure.figsize"] = (15,5)
        lab1 = []
        center_peak = []
        store_energy = []
        store_intensity = []
        final_filename_list = []
        filename_list = []
        for i in range(len(Total_length)):
            j = 0
            folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EELS/"
            directory = folder + Loc[i]
            data_folder = Path(directory)


            for filename in sorted(os.listdir(directory)):
                if filename.endswith(file) and 'thickness' not in filename and 'X' not in filename:
                    #print(filename)
                    f = open(data_folder / filename, 'r').read()
                    pure_list = f.split()
                    tot_intensity = [float(i) for i in pure_list]
                    spectrum_pure, background, spectrum_no_background = spectrum_EELS(tot_intensity)
                    print(i, j)
                    Energy_eV = np.linspace(0,Total_length[i][j], len(spectrum_pure)) #len of spectrum_EELS should be the same
                    Energy_eV_thickness = np.linspace(0,Total_length[i][j], len(tot_intensity)) #x-data for thickness measurments
                    norm = list(minmax_scale(np.array(spectrum_no_background))) #Choose which data_set one want, background, EELS
                    colors = ['darkred', 'goldenrod','lightskyblue', 'lightpink', 'black']
                    markers=[">", "^","s", "h", "+", "o", "v"]
                    final_filename = filename.replace("_", "\_")
                    print(final_filename[33:-4])

                    xGaus, norm_new, center, mu_1, sigma_1, y_peak, n_gaussians, rms, Gaussian_func, y_1 = gaussian_fit(Energy_eV,norm, num_gauss, range_index_gaus, peak, sigma)
                    dict_atoms = Atomic_percent('.txt', '/NEW_TEST.pdf')
                    n = 11

                    center_peak.append(center[sigma_1.index(max(sigma_1))])
                    store_energy.append(Energy_eV)
                    final_filename_list.append(final_filename[33:-4])
                    filename_list.append(filename[28:-4])
                    store_intensity.append(y_1)
                    '''
                    colors = np.flipud(plt.cm.copper(np.linspace(0,1,n)))#viridis

                    if final_filename[33:-4] in dict_atoms.keys():
                        print(dict_atoms[final_filename[33:-4]])
                        lab1.append('{:<22}:  $Ge_{{{:.2f}}}$$Se_{{{:.2f}}}$'.format(final_filename[33:-4], dict_atoms[final_filename[33:-4]][0], dict_atoms[final_filename[33:-4]][1]))
                    else:
                        lab1.append('{:<22}'.format(final_filename[33:-4]))



                    plt.plot(Energy_eV, y_1   ,color = colors[count])
                    plt.plot(Energy_eV, Gaussian_func[0])
                    plt.legend(lab1, loc='upper right', frameon=False)

                    plt.xlabel("Energy [eV]",fontsize=15)
                    plt.ylabel('Intensity [a.u]',fontsize=15)
                    plt.title('', fontsize=17)
                    plt.xlim([min(Energy_eV), max(Energy_eV)])
                    #plt.ylim([min(norm), max(norm) + max(norm)/10 ])
                    plt.tight_layout()
                    plt.savefig(folder + 'Combined' + ".pdf", format='pdf', dpi=1200)
                    '''
                    j += 1
                    count += 1


        sorted_final_filename_list = [x for _,x in sorted(zip(center_peak,final_filename_list))]
        sorted_filename_list = [x for _,x in sorted(zip(center_peak,filename_list))]

        sorted_store_energy = [x for _,x in sorted(zip(center_peak,store_energy))]
        sorted_store_intensity = [x for _,x in sorted(zip(center_peak,store_intensity))]

        sorted_center_peak = [x for x in sorted(center_peak)]

        for i in range(len(sorted_final_filename_list)):
            n = 11
            print(sorted_final_filename_list[i])
            print(sorted_center_peak[i])
            colors = np.flipud(plt.cm.copper(np.linspace(0,1,n)))#viridis
            if sorted_final_filename_list[i] in dict_atoms.keys():
                print(dict_atoms[sorted_final_filename_list[i]])
                lab1.append('{:<25}: $E_p$ {:.2f}:  $Ge_{{{:.2f}}}$$Se_{{{:.2f}}}$'.format(sorted_filename_list[i], sorted_center_peak[i], dict_atoms[sorted_final_filename_list[i]][0]/100, dict_atoms[sorted_final_filename_list[i]][1]/100))
            else:
                print(dict_atoms)
                lab1.append('{:<24} : $E_p$ {:.2f}:  $Ge_{{{:.2f}}}$$Se_{{{:.2f}}}$'.format(sorted_filename_list[i], sorted_center_peak[i], dict_atoms['GeSe\\_flake\\_ROI\\_1'][0]/100, dict_atoms['GeSe\\_flake\\_ROI\\_1'][1]/100))

            plt.plot(sorted_store_energy[i], sorted_store_intensity[i]   ,color = colors[i])
            plt.legend(lab1, loc='upper right', frameon=False, prop={'family': 'monospace', "size": 12})

            plt.xlabel("Energy [eV]",fontsize=15)
            plt.ylabel('Intensity [a.u]',fontsize=15)
            plt.title('', fontsize=17)
            plt.xlim([min(Energy_eV), max(Energy_eV)])
            #plt.ylim([min(norm), max(norm) + max(norm)/10 ])
            plt.tight_layout()
            plt.savefig(folder + 'Combined' + ".pdf", format='pdf', dpi=1200, bbox_inches='tight')


Total_length_GeSe_flake_height = [51.03, 252.87, 64.9, 179.82] #nm
Total_length_GeSe_flake_spectrum = [103.551, 103.03, 103.874, 103.504, 103.504] #nm

Total_length_GeSe_acu_spectrum = [103.905, 103.826, 103.807, 103.906] #nm

Total_length_GeSe_Needle_spectrum = [103.479, 103.628, 103.704, 103.828] #nm

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

sigma = [4.5,2.5,2.5, 2, 1, 1, 1] #Manually changing sigmas to fit to data...

folders = ['GeSe_flake','GeSe_accumulation','GeSe_Needle']
Total_length = [Total_length_GeSe_flake_spectrum, Total_length_GeSe_acu_spectrum,Total_length_GeSe_Needle_spectrum]

lines_peaks_plot = []
count_num_files = []
list_filenames = []

#Position('GeSe_Needle', '.txt',Total_length_GeSe_Needle_spectrum,  200,False, 'peak_table_plot',100, 1, sigma)
#for i in range(len(Total_length)):
#    EELS_Data_Spectrum(folders[i], '.txt',Total_length[i],  200,False, 'peak_table_plot',100, 1, sigma)

#peak_table_plots(folders[1], '.txt',Total_length[1],  200,False, 'peak_table_plot',100, 1, sigma)
total_EELS_data_spectrum(folders, '.txt',Total_length,  200,100, 1, sigma)
#EELS_Data_Spectrum(folders[0], '.txt',Total_length[0],  200,True, 'thickness',100, 1, sigma)
