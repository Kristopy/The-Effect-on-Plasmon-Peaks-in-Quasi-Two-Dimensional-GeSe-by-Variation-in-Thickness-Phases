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
from matplotlib import pyplot as plt
from astroML.datasets import fetch_vega_spectrum
from astroML.sum_of_norms import sum_of_norms, norm
from scipy.signal import find_peaks
import math
from cycler import cycler
import random

#rc('text', usetex=False)
#rc('font', **{'family': 'serif', 'serif': ['Random']})

def Reading_files(lines):
    #x = []
    #y = []

    a_dict = {}

    for i in range(len(lines)-1):
        if (i > 4 and i < 7):
            splitted_lines = lines[i].split()
            #x.append(splitted_lines[0])
            #y.append(float(splitted_lines[5]))  # second element u(x)
            a_dict[splitted_lines[0]] = [float(splitted_lines[1]), float(splitted_lines[2])]

    return a_dict

def Atomic_percent(file, save):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/EDX/"
    directory = folder
    data_folder = Path(directory)
    j = 0
    legends = []
    dict_atoms = {}
    plt.rcParams["figure.figsize"] = (15,5)

    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            final_filename = filename.replace("_", "\_")
            #print('-------------------------')
            #print('-------------------------')
            #print(final_filename[30:-4])
            f = open(data_folder / filename, 'r').read()
            at_dict = Reading_files(f.split("\n"))
            desired_order_list = ['Germanium', 'Selenium']
            reordered_dict = {k: at_dict[k] for k in desired_order_list}

            atoms1, at_percent1 = list(zip(*reordered_dict.items()))
            atoms, at_percent2 = list(atoms1), list(at_percent1)

            at_percent = [item[0] for item in at_percent2]
            Sigma = [item[1] for item in at_percent2]



            ratio_Se_Ge = at_percent[1]/at_percent[0]
            ratio_Ge_Se = at_percent[0]/at_percent[1]
            #print('Ratios: ', ratio_Ge_Se, ratio_Se_Ge)


            #print(final_filename[30:-4])
            dict_atoms[final_filename[30:-4]] = at_percent

            '''
            Weight  = [72.61, 78.96, 16.0] #Atomic weight of Ge, Se and O
            SUM = 0
            for i in range(len(Sigma)):

                SUM += Sigma[i]/Weight[i]

            Sigma_at = []
            for i in range(len(Sigma)):
                Sigma_at.append((Sigma[i]/Weight[i])/SUM)
            '''

            '''
            #print(atoms, '\n', at_percent, '\n', Sigma)

            x = np.linspace(0,30,2)


            #colors = ['dimgrey','darkred', 'goldenrod', 'orange','lightskyblue', 'lightpink', 'black', 'gold']
            n = 16
            colors = np.flipud(plt.cm.viridis(np.linspace(0,1,n)))#viridis
            markers=[">", "^","s", "h", "+", ".", "*", "1", "2", "3", "x", "d","D", "H", "8", "_"]
            if 'GeSe\_flake' in final_filename:
                print(j)
                plt.errorbar(x , at_percent, yerr=Sigma, linewidth = 0.7,linestyle=':',  color = colors[j],marker=markers[j], markersize=6)#label= '{:<16}-{:>8,f}'.format(final_filename[28:-4], at_percent[0]/(100), at_percent[1]/(100))) #$Ge_{%10.2f}Se_{%10.2f}$' % (final_filename[28:-4], at_percent[0]/(100), at_percent[1]/(100)))
            elif 'GeSe\_Acu' in final_filename:
                print(j)
                plt.errorbar(x , at_percent, yerr=Sigma, linewidth = 0.7,linestyle='dashdot',  color = colors[j],marker=markers[j], markersize=6)#label= '%s :        $Ge_{%10.2f}Se_{%10.2f}$' % (final_filename[28:-4], at_percent[0]/(100), at_percent[1]/(100)))
            else:
                print(j)
                plt.errorbar(x , at_percent, yerr=Sigma, linewidth = 0.7,linestyle='solid',  color = colors[j],marker=markers[j], markersize=6)#label= '%s :          $Ge_{%.2f}Se_{%.2f}$' % (final_filename[28:-4], at_percent[0]/(100), at_percent[1]/(100)))

                #:<28
            legends.append('{:<22}:  $Ge_{{{:.2f}}}$$Se_{{{:.2f}}}$'.format(filename[26:-4], at_percent[0]/(100), at_percent[1]/(100)))
            #plt.xlabel("",fontsize=15)
            legends_1 = [legends[j].replace(':   ', ':') for j in range(len(legends))]
            plt.ylabel('Atom concentration [At.$\\%$]',fontsize=15)
            plt.title("", fontsize=15)
            plt.legend(legends_1, loc='upper center', prop={'family': 'monospace'},ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.35) )
            plt.ylim([0, 100])
            plt.xticks([0,30],['Germanium', 'Selenium'])
            plt.savefig(directory + save, format='pdf', dpi=1200, bbox_inches='tight')

            j += 1
            '''
            'Comment out when not using, this py-file is being used in other files as well. Then one does not want to plot. '
            #return atoms, at_percent, Sigma, dict_atoms
    return dict_atoms

dict_atoms = Atomic_percent('.txt', '/NEW_TEST_1.pdf')
