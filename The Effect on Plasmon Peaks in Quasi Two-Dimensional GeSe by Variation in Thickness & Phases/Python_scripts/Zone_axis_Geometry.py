import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
import os
import itertools
from pathlib import Path
import sys
from matplotlib import pyplot as plt
import math
from side_by_side import print_side_by_side

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})

def Reading_files(lines):
    Zone_axis = []
    g_1_hkl = []
    g_2_hkl = []

    Ratio_g_1_g_2 = []
    Ratio_g_2_g_1 = []

    angle = []
    for i in range(len(lines)-1):
        if (i > 1 and i != 0 ):
            splitted_lines = lines[i].split()

            Zone_axis.append(' '.join(splitted_lines[0:3]))  # First element x-values corrensponding to Zone Axis
            g_1_hkl.append(' '.join(splitted_lines[3:6]))  # hkl of g_1 for each zone axis
            g_2_hkl.append(' '.join(splitted_lines[8:11]))  # hkl of g_2 for each zone axis
            Ratio_g_1_g_2.append(float(splitted_lines[13]))  # second e
            Ratio_g_2_g_1.append(float(splitted_lines[14]))  # second e
            angle.append(float(splitted_lines[15]))  # second e

    return Zone_axis, g_1_hkl, g_2_hkl, Ratio_g_1_g_2, Ratio_g_2_g_1, angle

def Position_1(file, dir, b_1, b_2, b_3, error):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/ZoneAxis_Diffraction_pattern/" + dir
    #sys.stdout = open("/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/Txt_file_ZoneAxis/2" + file, "w")

    directory = folder
    data_folder = Path(directory)

    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):

            print('-------------------------------------------')
            print('-------------------------------------------')
            print('FILENAME: ', filename)
            print('-------------------------------------------')

            f = open(data_folder / filename, 'r').read()
            Zone_axis, g_1_hkl, g_2_hkl, Ratio_g_1_g_2, Ratio_g_2_g_1, angle = Reading_files(f.split("\n"))


            #Calculation different ration from high resolution STEM images. Three different vectors, b1 b2 b3.
            # Ratio of one of these should correlate with some of the Zone axises.
            Ratio_b_1_b_2 = b_1/b_2
            Ratio_b_2_b_1 = b_2/b_1

            Ratio_b1_b_3 = b_1/b_3
            Ratio_b3_b_1 = b_3/b_1

            Ratio_b2_b_3 = b_2/b_3
            Ratio_b3_b_2 = b_3/b_2
            Ratio_b3_b_3 = b_3/b_3
            for i in range(len(Zone_axis)):
                if angle[i] == 90 and abs(Ratio_g_1_g_2[i]-Ratio_b_1_b_2) <= 0.1:
                    print(Zone_axis[i], ':', Ratio_g_1_g_2[i], Ratio_g_2_g_1[i], angle[i], 'b1/b2, b2/b1: ', Ratio_b_1_b_2, Ratio_b_2_b_1)


                if angle[i] == 90 and abs(Ratio_g_1_g_2[i]-Ratio_b3_b_3) <= 0.1:
                    print(Zone_axis[i], ':', Ratio_g_1_g_2[i], Ratio_g_2_g_1[i], angle[i], 'b3/b3: ', Ratio_b3_b_3)
                '''
                if abs(Ratio_g_1_g_2[i] - Ratio_b_1_b_2) <= error or abs(Ratio_g_2_g_1[i] - Ratio_b_2_b_1) <= error:
                    print('---------------------------------------------')
                    print('Zone Axis for b1/b2 and b2/b1 : ' f'{Zone_axis[i]} \n\n',
                    'hkl plane for g1: ', f'{g_1_hkl[i]} \n',
                    'hkl plane for g2: ', f'{g_2_hkl[i]} \n \n',
                    'Ratio: g1/g2: ', f'{Ratio_g_1_g_2[i]}', '\n        g2_g1: ',f'{Ratio_g_2_g_1[i]}\n \n',
                    'Ratio: b1/b2: ', f'{Ratio_b_1_b_2}', '\n        b2_b1: ',f'{Ratio_b_2_b_1}\n\n',
                    'Difference: b1/b2: ', f'{abs(Ratio_g_1_g_2[i] - Ratio_b_1_b_2)}', '\n             b2_b1: ',f'{abs(Ratio_g_2_g_1[i] - Ratio_b_2_b_1) }\n')

                if abs(Ratio_g_1_g_2[i] - Ratio_b1_b_3) <= error or abs(Ratio_g_2_g_1[i] - Ratio_b3_b_1) <= error:
                    print('---------------------------------------------')
                    print('Zone Axis for b1/b3 and b3/b1 : ' f'{Zone_axis[i]} \n\n' ,
                    'hkl plane for g1: ', f'{g_1_hkl[i]} \n',
                    'hkl plane for g2: ', f'{g_2_hkl[i]} \n \n',
                    'Ratio: g1/g2: ', f'{Ratio_g_1_g_2[i]}', '\n        g2_g1: ',f'{Ratio_g_2_g_1[i]}\n\n',
                    'Ratio: b1/b3: ', f'{Ratio_b1_b_3}', '\n        b3_b1: ',f'{Ratio_b3_b_1}\n\n',
                    'Difference: b1/b3: ', f'{abs(Ratio_g_1_g_2[i] - Ratio_b1_b_3)}', '\n             b3_b1: ',f'{abs(Ratio_g_2_g_1[i] - Ratio_b3_b_1)}\n')


                if abs(Ratio_g_1_g_2[i] - Ratio_b2_b_3) <= error or abs(Ratio_g_2_g_1[i] - Ratio_b3_b_2) <= error:
                    print('---------------------------------------------')
                    print('Zone Axis for b2/b3 and b3/b2 : ' f'{Zone_axis[i]} \n\n' ,
                    'hkl plane for g1: ', f'{g_1_hkl[i]} \n',
                    'hkl plane for g2: ', f'{g_2_hkl[i]} \n \n',
                    'Ratio: g1/g2: ', f'{Ratio_g_1_g_2[i]}', '\n        g2_g1: ',f'{Ratio_g_2_g_1[i]}\n \n',
                    'Ratio: b2/b3: ', f'{Ratio_b2_b_3}', '\n        b3_b2: ',f'{Ratio_b3_b_2}\n\n',
                    'Difference: b2/b3: ', f'{abs(Ratio_g_1_g_2[i] - Ratio_b2_b_3)}', '\n             b3_b2: ',f'{abs(Ratio_g_2_g_1[i] - Ratio_b3_b_2)}\n')
                    '''

b_1 = 5.235#5.23560
b_2 = 4.62 #4.73933
b_3 = 3.52#7.05#3.52113


files = ['Pnma.txt', 'Pbma.txt', 'Pcma.txt', '2.txt']
loc = ['Pnma', 'Pbma', 'Pcma', 'GeSe2']

for i in range(len(files)):
    Position_1('Zone_axis_geometry_GeSe_' + files[i], loc[i], b_1, b_2, b_3,0.1     )

#sys.stdout.close()
