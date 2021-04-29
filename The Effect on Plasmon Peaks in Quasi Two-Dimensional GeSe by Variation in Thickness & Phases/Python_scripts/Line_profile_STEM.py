import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
import os
import itertools
from pathlib import Path
import sys
from matplotlib import pyplot as plt
import math
from scipy.signal import find_peaks
from sklearn.preprocessing import minmax_scale
from astroML.sum_of_norms import sum_of_norms, norm
from scipy.signal import find_peaks, peak_widths
from EELS_Data_Spectrum import gaussian_fit
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})




def Position_1(file, Total_length, Number_of_peaks, Lenght_peaks):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/STEM/GeSe_flake/GeSe_STEM/"
    directory = folder
    data_folder = Path(directory)
    i = 0
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            print('-------------------------------------------')
            print('-------------------------------------------')
            print(filename)
            print('-------------------------------------------')
            print('-------------------------------------------')
            f = open(data_folder / filename, 'r').read()
            pure_list = f.split()

            intensity = [float(i) for i in pure_list]
            norm = list(minmax_scale(np.array(intensity)))
            print(norm)
            Lenght_x = np.linspace(0,Total_length, len(intensity))

            l_b_p = Lenght_peaks/(Number_of_peaks-1)
            colors = ['darkred', 'goldenrod','lightskyblue', 'lightpink', 'black']
            plt.rcParams["figure.figsize"] = (15,10)
            plt.figure()

            plt.fill_between(Lenght_x[5:165], norm[5:165] , alpha=0.4, color = colors[i], label= 'Intensity STEM, num: %s' % file[1:2])

            plt.legend(loc='upper right', prop={"size": 13}, frameon=False)

            textstr = '\n'.join((
            r'Total lenght = %.3fnm' % (Total_length, ),
            r'Number of peaks= %.0f' % (Number_of_peaks, ),
            r'Lenght between peak 1 - %.0f : %.3fnm' % (Number_of_peaks, Lenght_peaks),
            r'Lenght between peaks = %.3fnm' % (l_b_p, )))

            #props = dict(boxstyle='italic',facecolor = 'grey', alpha =  0.2, pad = 10)
            props = dict({'facecolor': 'grey', 'alpha': 0.2, 'pad': 10})
            plt.text(max(Lenght_x) - max(Lenght_x)/4.6, max(intensity) + max(intensity)/20, textstr, fontsize=12,
            verticalalignment='top', bbox=props)

            plt.xlabel("Lenght [nm]",fontsize=15)
            plt.ylabel('Intensity [a.u]',fontsize=15)
            plt.title('', fontsize=17)
            plt.xlim([min(Lenght_x[5:]), max(Lenght_x[:165])])
            plt.ylim([min(norm), max(norm) + max(norm)/10 ])
            plt.tight_layout()

            plt.savefig(folder + file[0:2] + "test.pdf", format='pdf', dpi=1200)

            i += 0

files = ['#1.txt', '#2.txt', '#3.txt',  '#4.txt']
Number_of_peaks = [19, 17, 11, 9]
Total_length = [5.4205, 4.9544, 4.2014, 3.903] #nm
Lenght_peaks = [5.1048, 4.6294, 3.8155, 3.4294] #nm

Position_1('#1.txt', Total_length[0], Number_of_peaks[0], Lenght_peaks[0])
#for i in range(len(files)):
#    Position_1(files[i], Total_length[i], Number_of_peaks[i], Lenght_peaks[i], i)
