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
import pandas as pd

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})

def Reading_files(lines):
    Zone_axis = []
    g_1_hkl = []
    g_2_hkl = []

    Ratio_g_1_g_2 = []
    Ratio_g_2_g_1 = []
    for i in range(len(lines)-1):
        if (i > 1 and i != 0 ):
            splitted_lines = lines[i].split()

            Zone_axis.append(' '.join(splitted_lines[0:3]))  # First element x-values corrensponding to Zone Axis
            g_1_hkl.append(' '.join(splitted_lines[3:6]))  # hkl of g_1 for each zone axis
            g_2_hkl.append(' '.join(splitted_lines[8:11]))  # hkl of g_2 for each zone axis
            Ratio_g_1_g_2.append(float(splitted_lines[13]))  # second e
            Ratio_g_2_g_1.append(float(splitted_lines[14]))  # second e

    return Zone_axis, g_1_hkl, g_2_hkl, Ratio_g_1_g_2, Ratio_g_2_g_1


def Position_1(file, dir, b_1, b_2, b_3, error):
    folder="/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/ZoneAxis_Diffraction_pattern/" + dir
    #sys.stdout = open("/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/Txt_file_ZoneAxis/" + file, "w")

    directory = folder
    data_folder = Path(directory)

    for filename in sorted(os.listdir(directory)):
        if filename.endswith(file):
            print('-------------------------------------------')
            print('-------------------------------------------')
            print('FILENAME: ', filename)
            print('-------------------------------------------')
            list_filename.append(filename)
            #print(list_filename)
            f = open(data_folder / filename, 'r').read()
            Zone_axis, g_1_hkl, g_2_hkl, Ratio_g_1_g_2, Ratio_g_2_g_1 = Reading_files(f.split("\n"))

            #Calculation different ration from high resolution STEM images. Three different vectors, b1 b2 b3.
            # Ratio of one of these should correlate with some of the Zone axises.
            Ratio_b_1_b_2 = b_1/b_2
            Ratio_b_2_b_1 = b_2/b_1

            Ratio_b1_b_3 = b_1/b_3
            Ratio_b3_b_1 = b_3/b_1

            Ratio_b2_b_3 = b_2/b_3
            Ratio_b3_b_2 = b_3/b_2



            for i in range(len(Zone_axis)):
                if abs(Ratio_g_1_g_2[i] - Ratio_b_1_b_2) <= error or abs(Ratio_g_2_g_1[i] - Ratio_b_2_b_1) <= error:
                    print('---------------------------------------------')
                    print('Zone Axis for b1/b2 and b2/b1 : ' f'{Zone_axis[i]} \n\n' ,
                    'hkl plane for g1: ', f'{g_1_hkl[i]} \n',
                    'hkl plane for g2: ', f'{g_2_hkl[i]} \n \n',
                    'Ratio: g1/g2: ', f'{Ratio_g_1_g_2[i]}', '\n        g2_g1: ',f'{Ratio_g_2_g_1[i]}\n \n',
                    'Ratio: b1/b2: ', f'{Ratio_b_1_b_2}', '\n        b2_b1: ',f'{Ratio_b_2_b_1}\n')

                if abs(Ratio_g_1_g_2[i] - Ratio_b1_b_3) <= error or abs(Ratio_g_2_g_1[i] - Ratio_b3_b_1) <= error:
                    print('---------------------------------------------')
                    print('Zone Axis for b1/b3 and b3/b1 : ' f'{Zone_axis[i]} \n\n' ,
                    'hkl plane for g1: ', f'{g_1_hkl[i]} \n',
                    'hkl plane for g2: ', f'{g_2_hkl[i]} \n \n',
                    'Ratio: g1/g2: ', f'{Ratio_g_1_g_2[i]}', '\n        g2_g1: ',f'{Ratio_g_2_g_1[i]}\n \n',
                    'Ratio: b1/b3: ', f'{Ratio_b1_b_3}', '\n        b3_b1: ',f'{Ratio_b3_b_1}\n')


                if abs(Ratio_g_1_g_2[i] - Ratio_b2_b_3) <= error or abs(Ratio_g_2_g_1[i] - Ratio_b3_b_2) <= error:
                    print('---------------------------------------------')
                    print('Zone Axis for b2/b3 and b3/b2 : ' f'{Zone_axis[i]} \n\n' ,
                    'hkl plane for g1: ', f'{g_1_hkl[i]} \n',
                    'hkl plane for g2: ', f'{g_2_hkl[i]} \n \n',
                    'Ratio: g1/g2: ', f'{Ratio_g_1_g_2[i]}', '\n        g2_g1: ',f'{Ratio_g_2_g_1[i]}\n \n',
                    'Ratio: b2/b3: ', f'{Ratio_b2_b_3}', '\n        b3_b2: ',f'{Ratio_b3_b_2}\n')
    return list_filename



b_1 = 5.23560
b_2 = 4.73933
b_3 = 3.52113


files = ['Pnma.txt', 'Pbma.txt', 'Pcma.txt', '2.txt']
loc = ['Pnma', 'Pbma', 'Pcma', 'GeSe2']
list_filename = []
Headers = ['Zone Axis', 'hkl g1', 'hkl g2', '[g1/g2]', '[g2/g1]']

for i in range(len(files)):
    list_filename = Position_1('Zone_axis_geometry_GeSe_' + files[i], loc[i], b_1, b_2, b_3, 0.003)




#list_filename =  [ '','','Zone_axis_geometry_GeSe_Pnma.txt','','','' ,'','Zone_axis_geometry_GeSe_Pbma.txt','','','','', 'Zone_axis_geometry_GeSe_Pcma.txt','','','','', 'Zone_axis_geometry_GeSe_2.txt','','']

list_filename =  ['Zone_axis_geometry_GeSe_Pnma.txt','Zone_axis_geometry_GeSe_Pbma.txt','Zone_axis_geometry_GeSe_Pcma.txt','Zone_axis_geometry_GeSe_2.txt']

print(len(list_filename))
pd.set_option("max_colwidth", 10)
pd.set_option("colheader_justify", "left")

df2 = pd.DataFrame(np.array([Headers*4, Headers*4, Headers*4]),
                   columns=list_filename)

df2.to_csv('/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/TEM/ZoneAxis_Diffraction_pattern/OK.csv')
#sys.stdout.close()
'''

df= pd.DataFrame({'Table of 9': [9,18,27],
        'Table of 10': [10,20,30]})

for i in range(1,11):
    df=df.append({'Table of 9':i*9,'Table of 10':i*10},ignore_index=True)

print(df)
'''
