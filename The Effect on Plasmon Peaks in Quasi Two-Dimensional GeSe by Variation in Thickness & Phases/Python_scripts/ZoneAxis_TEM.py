import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
import os
import itertools
from pathlib import Path
import sys
from matplotlib import pyplot as plt
import math


rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})

def Reading_files(lines):
    x = []
    y_float = []

    for i in range(len(lines)-1):
        if (i > 4 and i != 0 ):
            splitted_lines = lines[i].split(',')
            x.append(' '.join(splitted_lines[0:3]))  # First element x-values corrensponding to miller indices plane.
            y_float.append(splitted_lines[8])  # second e
    y  = [float(i) for i in y_float]

    return x,y

def Position_1(file, dir, b_1, b_2, b_3):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/ZoneAxis_Diffraction_pattern/" + dir
    #sys.stdout = open("/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/Txt_file_ZoneAxis/" + file, "w")


    directory = folder
    data_folder = Path(directory)

    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            print('-------------------------------------------')
            print('-------------------------------------------')
            print(filename)
            print('-------------------------------------------')
            print('-------------------------------------------')
            f = open(data_folder / filename, 'r').read()
            indices,abs_g = Reading_files(f.split("\n"))


            min_r = b_1/b_2
            max_r = b_1/b_3

            #1.34375, 1.492146
            for i in range(25):
                if abs_g[i] != abs_g[i+1]:
                    print(f'{indices[i]:15} ==> {abs_g[i]:10f} \n' )
                    ratio = abs_g[0]/abs_g[i+1]
                    print('Ratio: ',f'{indices[0]}/{indices[i+1]} :', ratio)
                    if abs(ratio - min_r) <= 0.14 or abs(ratio - max_r) <= 0.15:
                        print('-------------------------------------------')
                        print('-------------------------------------------')
                        print('************EUREKAAAA :',f'{i}', f'{indices[0]}/{indices[i+1]}************ ',abs(ratio - min_r), abs(ratio - max_r)  )
                        print('-------------------------------------------')
                        print('-------------------------------------------')
                else:
                    print(f'{indices[i]:15} ==> {abs_g[i]:10f}' )


            print('-------------------------------------------')

            print('Comparing Diffraction pattern miller indices to real STEM FFT values:')


            for i in range(8):
                if abs(abs_g[i] - b_1) <= 0.01:
                    print(f'{indices[i]:15} ==> {abs_g[i]:10f}' )
                elif abs(abs_g[i] - b_2) <= 0.01:
                    print(f'{indices[i]:15} ==> {abs_g[i]:10f}' )
                elif abs(abs_g[i] - b_3) <= 0.01:
                    print(f'{indices[i]:15} ==> {abs_g[i]:10f}' )



files = ['Pnma.txt', 'Pbma.txt', 'Pcma.txt']
loc = ['Pnma', 'Pbma', 'Pcma', 'GeSe2']

Position_1('Zone_axis_geometry_GeSe_Pnma.txt', 'Pnma', 0.285, 0.211, 0.191)
#sys.stdout.close()
