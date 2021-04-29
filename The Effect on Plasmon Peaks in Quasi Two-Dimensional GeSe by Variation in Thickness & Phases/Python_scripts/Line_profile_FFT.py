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

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})




def Position_1(file, Total_length, i, dir):
    #folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/STEM/GeSe_flake/GeSe_STEM/"
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/STEM/GeSe_flake/FFT/" + dir
    directory = folder
    data_folder = Path(directory)
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            print('-------------------------------------------')
            print('-------------------------------------------')
            final_filename = filename.replace("_", " ")
            label_add = final_filename.replace(" ", '-' )
            print(filename)
            print('-------------------------------------------')
            print('-------------------------------------------')
            f = open(data_folder / filename, 'r').read()
            pure_list = f.split()
            labels =['$\\bf{b_1}$', '$\\bf{b_2}$', '$\\bf{b_3}^1$', '$\\bf{b_3}^2 $']
            intensity = [abs(float(i)) for i in pure_list]

            Lenght_x = np.linspace(0,Total_length, len(intensity))
            print(intensity[8])
            colors = ['darkred', 'goldenrod','lightskyblue', 'lightpink', 'black']

            plt.rcParams["figure.figsize"] = (15,5)
            plt.figure()
            print(i)
            plt.yscale('log')
            plt.fill_between(Lenght_x, intensity , alpha=0.4, color = colors[i], label= r'Intensity FFT, %s' %labels[i] )
            #plt.plot(Lenght_x, intensity, color = colors[i], label= 'Intensity FFT, num')

            plt.legend(loc='upper right', prop={"size": 14}, frameon=False)

            textstr = '\n'.join((
            r'Total lenght = %.3fnm' % (Total_length, )))

            #props = dict(boxstyle='italic',facecolor = 'grey', alpha =  0.2, pad = 10)
            props = dict({'facecolor': 'grey', 'alpha': 0.2, 'pad': 10})
            plt.text(max(Lenght_x) , max(intensity) , textstr, fontsize=12,
            verticalalignment='top', bbox=props)

            plt.xlabel("Lenght [1/nm]",fontsize=15)
            plt.ylabel('Intensity [a.u]',fontsize=15)
            plt.title("", fontsize=17)
            plt.xlim([min(Lenght_x), max(Lenght_x)])
            plt.ylim([min(intensity), max(intensity) + max(intensity)/10 ])
            plt.tight_layout()

            plt.savefig(folder + filename[13:-4] + ".png", format='png', dpi=1200)


'''
#For basic line profile of high resolution image STEM
files = ['#1.txt', '#2.txt', '#3.txt']
Number_of_peaks = [19, 17, 11]
Total_length = [5.4205, 4.9544, 4.2014] #nm
Lenght_peaks = [5.1048, 4.6294, 3.8155] #nm
'''

#Line profile of FFT bot 15.10.03 and main high resoulution image 15.21.52

loc = ['15.10.03', '15.21.52']

files = ['b_1.txt', 'b_2.txt', 'b_3_1.txt', 'b_3_2.txt']

Total_length_1 = [11.58,19.6, 14.93, 15.75] #nm
peak_lenght_1 = []
Total_length_2 = [12.03, 11.58, 14.72,15.37] #nm


for i in range(len(files)):
    Position_1(files[i], Total_length_1[i], i, '15.10.03')
    #Position_1(files[i], Total_length_2[i], i, '15.21.52')
