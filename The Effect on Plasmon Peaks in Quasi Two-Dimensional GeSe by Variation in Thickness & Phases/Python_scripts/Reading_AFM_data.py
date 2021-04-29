import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib import rc, rcParams
from scipy.stats import norm
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Random']})

directory = "/Users/kristoffervarslott/Documents/MENA/Master_MENA/Masteroppgave/Experimental_data/"

path_file_1 = directory + "AFM-Analysis/2D-GeSe_#4"
path_file_2 = directory + "AFM-Analysis/2D-GeSe_#4"
path_file_3 = directory + "Optical_Microscope/GeSe_Scale"
path_file_4 = directory + "Optical_Microscope/GeSe_Scale/ImageJ_100X"

graph_data_1 = open(path_file_1 + '/GeSe_#4_006_Texfile_HeightLine.txt', 'r').read()
graph_data_2 = open(path_file_2 + '/GeSe_#4_pos_2_002_Texfile_HeightLine.txt', 'r').read()
graph_data_3 = open(path_file_3 + '/Distribution_flakes_GeSe_#4_Si_substrate_tot.txt', 'r').read()
graph_data_4 = open(path_file_4 + '/Contrast_value2.txt', 'r').read()

#Simple reading files function - be carefull, dimension may vary in original data, and this
#needs to be handled in this file - also, there are 6 variables created along the x-coordinated.
#Reads 1 - lines-1 and six colums.
def Reading_files(lines):
    x = []
    y = []

    x1 = []
    y1 = []

    x2 = []
    y2 = []

    for i in range(len(lines)-1):
        if (i >= 3 and i < len(lines)-1):
            splitted_lines = lines[i].split()

            if '-' not in splitted_lines:
                c = list(map(float, splitted_lines))
                x.append(c[0]*10**6)  # First element x-values
                y.append(c[1]*10**9)  # second element u(x)
                if len(c) > 2:
                    x1.append(c[2]*10**6)  # First element x-values
                    y1.append(c[3]*10**9)  # second element u(x)
                    x2.append(c[4]*10**6)  # First element x-values
                    y2.append(c[5]*10**9)  # second element u(x)
    return x,y,x1,y1,x2,y2


#Reading and plotting AFM data - line profiles extrated from raw-file.
#Two positions defined, Pos_1 for simple 1 figure plots where the lines coincides
#Pos_2 is used when the scale is off between the line profiles, two subplots are created here.
def Position_1():
    fig, ax = plt.subplots()

    #Calculating the avergae increase/decline in thickness.

    #Line 1) one step.
        # top goes from x[0:26] find average y value.
        # bottom goes from x[31:] find average y value.

    line_1_y_mean_1 = np.mean(y[0:26])
    line_1_y_mean_2 = np.mean(y[31:])

    #line 2) one step.
        # top goes from x[0:26] find average y value.
        # bottom goes from x[28:] find average y value.

    line_2_y_mean_1 = np.mean(y1[0:26])
    line_2_y_mean_2 = np.mean(y1[28:])

    #line 3) three steps. Four lines from top to bottom.

        #  x[0:14] find average y value.
        #  x[19:27] find average y value.

        #  x[30:36] find average y value.
        #  x[41:]  find average y value.

    line_3_y_mean_1 = np.mean(y2[0:14])
    line_3_y_mean_2 = np.mean(y2[19:27])
    line_3_y_mean_3 = np.mean(y2[30:36])
    line_3_y_mean_4 = np.mean(y2[41:])

    height_1 =[]
    height_2=[]
    height_3=[]
    height_list_1= [line_1_y_mean_1,line_1_y_mean_2]
    height_list_2= [line_2_y_mean_1,line_2_y_mean_2]
    height_list_3= [line_3_y_mean_1,line_3_y_mean_2,line_3_y_mean_3,line_3_y_mean_4]
    for i in range(len(height_list_1)-1):
        height_1.append(height_list_1[i] - height_list_1[i+1])
    for i in range(len(height_list_2)-1):
        height_2.append(height_list_2[i] - height_list_2[i+1])
    for i in range(len(height_list_3)-1):
        height_3.append(height_list_3[i] - height_list_3[i+1])


    #print(height_1,height_2,height_3)


    x_values = [x[26],x[31], x1[26], x1[28],x2[19], x2[30], x2[41]]
    #print(x_values)
    l = 0
    for i in range(3):
        tot_height = [[line_1_y_mean_1,line_1_y_mean_2], [line_2_y_mean_1,line_2_y_mean_2], [line_3_y_mean_1,line_3_y_mean_2,line_3_y_mean_3,line_3_y_mean_4]]
        #print(tot_height[i])

        for j in range(len(tot_height[i])-1):

            mean_value = (tot_height[i][j+1]-tot_height[i][j])/2 + tot_height[i][j]
            #print(mean_value)
            print(tot_height[i][j]-tot_height[i][j+1])
            ax.plot([x_values[l],x_values[l]],[tot_height[i][j],tot_height[i][j+1]], linewidth=0.5)

            #print(tot_height[i][j]-tot_height[i][j+1])
            #print((tot_height[i][j]-tot_height[i][j+1])/2 + tot_height[i][j])
            ax.annotate('%.2fnm' % (tot_height[i][j]-tot_height[i][j+1]),
                        xy=(x_values[l],mean_value),
                        xytext=(15, 0),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')


            l +=1

    plt.plot(x, y, linewidth=0.7,color= 'dimgrey',marker= "o", markersize=3, label="Line 1")
    plt.plot(x1, y1, linewidth=0.7,color= 'darkred',marker= "^", markersize=3, label="Line 2")
    ax.plot(x2, y2, linewidth=0.7,color= 'goldenrod',marker=  "s",markersize=3, label="Line 3")

    ax.legend(loc='upper right', prop={"size": 13}, frameon=False)


    ax.set_xlabel("Position along x,y-plane [$\\mu$m]",fontsize=15)
    ax.set_ylabel('Height [nm]',fontsize=15)

    plt.title("Height profile of GeSe; Position I", fontsize=17)

    ax.set_xlim(0,0.5)

    plt.savefig(path_file_1 + "/TEST_height.eps", format='eps', dpi=1200)

def Position_2():
    fig, axis = plt.subplots(3,1)


    line_1_y_mean_1 = np.mean(y[0:26])
    line_1_y_mean_2 = np.mean(y[30:])

    print(line_1_y_mean_1-line_1_y_mean_2)


    line_2_y_mean_1 = np.mean(y1[0:21])
    line_2_y_mean_2 = np.mean(y1[24:])

    print(line_2_y_mean_1-line_2_y_mean_2)


    line_3_y_mean_1 = np.mean(y2[0:26])
    line_3_y_mean_2 = np.mean(y2[29:])

    print(line_3_y_mean_1-line_3_y_mean_2)

    axis[0].plot([x[26],x[26]],[line_1_y_mean_1,line_1_y_mean_2], linewidth=0.5)

    axis[1].plot([x1[21],x1[21]],[line_2_y_mean_1,line_2_y_mean_2], linewidth=0.5)

    axis[2].plot([x2[25],x2[25]],[line_3_y_mean_1,line_3_y_mean_2], linewidth=0.5)


    axis[0].plot([x[30],x[30]],[line_1_y_mean_1,line_1_y_mean_2], linewidth=0.5)

    axis[1].plot([x1[24],x1[24]],[line_2_y_mean_1,line_2_y_mean_2], linewidth=0.5)

    axis[2].plot([x2[29],x2[29]],[line_3_y_mean_1,line_3_y_mean_2], linewidth=0.5)

    axis[0].plot(x, y, linewidth=0.7,color= 'dimgrey',marker= "o", markersize=3, label="Line 1")
    axis[1].plot(x1, y1, linewidth=0.7,color= 'darkred',marker= "^", markersize=3, label="Line 2")
    axis[2].plot(x2, y2, linewidth=0.7,color= 'goldenrod',marker=  "s",markersize=3, label="Line 3")

    axis[0].legend(loc='upper right', prop={"size": 13}, frameon=False)
    axis[1].legend(loc='upper right', prop={"size": 13}, frameon=False)
    axis[2].legend(loc='upper right', prop={"size": 13}, frameon=False)

    axis[2].set_xlabel("Position along x,y-plane [$\\mu$m]",fontsize=15)
    axis[1].set_ylabel('Height [nm]',fontsize=15)

    axis[0].set_title("Height profile of GeSe; Position II", fontsize=17)

    axis[0].get_xaxis().set_visible(False)
    axis[1].get_xaxis().set_visible(False)

    axis[0].set_yticks(np.arange(min(y), max(y), 1))
    axis[1].set_yticks(np.arange(min(y1), max(y1), 1))
    axis[2].set_yticks(np.arange(min(y2), max(y2), 1))

    axis[0].set_xlim(0,0.4)
    axis[1].set_xlim(0,0.4)
    axis[2].set_xlim(0,0.4)

    #axis[2].autoscale(enable=True, axis='x', tight=True)
    plt.savefig(path_file_2 + "/TEST.eps", format='eps', dpi=1200)

#Normal distribution of flakes - Showing the average size of flakes from GeSe dispensed in isopropanol.
#Potentially this could be done for 2D+flakes only.
def Normal_distribution():
    plt.rcParams["figure.figsize"] = (10,5)
    plt.figure()
    mu, sigma = abs(np.mean(y1)),abs(np.std(y1, ddof=1))
    #mu, sigma = norm.fit(y1)

    #bins:34
    #weights=np.ones(len(y1)) / len(y1), implement in plt.hist if sum of bins should be 1 and set
    count, bins, ignored = plt.hist(y1,  bins=42, density=True,color='darkred',  alpha=0.7, label='Histogram')

    xmin, xmax = plt.xlim()
    xGaus = np.linspace(xmin, xmax, 100)
    Gaus = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (xGaus - mu)**2 / (2 * sigma**2))
    plt.plot(xGaus, Gaus,
         linewidth=1, linestyle='--', color='goldenrod', label='Gaussian distribution')

    #xmin, xmax = plt.xlim()
    #x = np.linspace(xmin, xmax, 100)
    #p = norm.pdf(x, mu, sigma)
    #plt.plot(x, p, 'k', linewidth=2)

    print('Mean=%.3f, Standard Deviation=%.3f' % (mu, sigma))
    print(sorted(y1))


    plt.legend(loc='upper right', prop={"size": 13}, frameon=False)
    plt.xlabel("Lenght along x,y-plane [$\\mu$m]",fontsize=15)
    plt.ylabel('Density',fontsize=15)
    #plt.title("Height profile of GeSe; Position I", fontsize=17)
    plt.text(mu*3.119 , max(count)*0.7, 'Mean=%.3f \n Standard Deviation=%.3f' % (mu, sigma), style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.2, 'pad': 10}, multialignment="left")
    plt.xticks(np.arange(0, max(xGaus), 1.0))

    plt.xlim(0,max(xGaus))
    plt.tight_layout()
    plt.savefig(path_file_3 + "/histogram.pdf", format='pdf', dpi=1200)


#Working progress: Determining thickness through contrast anaysis.
def Contrast_analysis():
    plt.figure()

    plt.plot(x, y, linewidth=0.7,color= 'dimgrey',marker= "o", markersize=3, label="Line 1")

    #print(np.mean(y))

    #Mean value of substrate contrast = 157.59061297935102

    plt.plot([0,max(x)],[max(y),max(y)], linewidth=0.7, linestyle='--', color='r')
    plt.plot([0,max(x)],[min(y),min(y)], linewidth=0.7, linestyle='--', color='r')
    plt.legend(loc='upper right', prop={"size": 13}, frameon=False)
    plt.xlabel("Position along x,y-plane [$\\mu$m]",fontsize=15)
    plt.ylabel('Contrast',fontsize=15)
    plt.title("Height profile of GeSe; Position I", fontsize=17)
    plt.savefig(path_file_4 + "/Contrast_plot.eps", format='eps', dpi=1200)


x,y,x1,y1,x2,y2 = Reading_files(graph_data_2.split("\n"))
#Position_1()
Position_2()
#Normal_distribution()
#Contrast_analysis()
