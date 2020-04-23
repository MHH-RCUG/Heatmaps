
#!/usr/bin/env python
# coding: utf-8


#Author: Konstantinos Sifakis
#Script made: Jan 2020-Feb 2020
#Script made for: MHH - RCUG


#________CHECK what is needed:
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import os
from pylab import savefig
import seaborn as sns; sns.set(color_codes=True)
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, linkage
#import warnings
from scipy.cluster.hierarchy import dendrogram
from argparse import RawTextHelpFormatter   #used in Description in order to help with text format.



#warnings.filterwarnings("ignore",category=UserWarning)
#Needed to read file  [ --filename ]  and then [ file directory ] so that file would open as a tsv file

parser = argparse.ArgumentParser(description="""   Create heatmaps/correlation heatmaps and dendrograms generated from csv
files that are generated from wochenende pipeline. Those files should be filtered, sorted and be into a csv format as the above.

Example: MCF01_S7_R1.ndp.lc.trm.s.mq30.01mm.dup.calmd.bam.txt.filt.sort.csv

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Before running python script, make sure that you have a .txt file that incudes all the names-samples that you wish to combine, in order to create heatmaps.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
This .txt file should be used an argument --batch_files .txt (location)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   Each run of this script creates several heatmaps and dendrograms. In order to prevent overwriting, there are 2 different implemented label sets per .png files.

   1) A number-label that the user changes in every run that he/she decides.(input)
   2) A set of numbers that the original input for each argument user used or Default values of the missing arguments, so that names of each heatmap can give information for that heatmap.

    Each file has read counts into 1 column and those columns are merged into a matrix.
For instance, a matrix has rows that are bacterial strains and columns that are files-sample names.
The ammount of data that are used are divided into Quarters for rows and columns too!
        User can choose how many Quarters of the files-samples COLUMNS should be used for heatmap. Subsets can be set from argument (( filter_sStart )) and values according to limit are:
minimum:0
25%:1,
50%:2
75%:3
maximum:4

        User can choose any subset with those limits as starting poing and as end point (( filter_sEnd ))
Of course it is suggested to have (( filter_sEnd )) > (( filter_sStart )) and (( filter_sEnd )) should be always set at value:4 in order to have the most informative files-samples

    Another filtering is done automatically and creates multiple heatmaps depending on the ROWS filtering per quartiles this time.
Different dimension on matrix for those Quarters then before, this time filtering is on ROWS and no input needed.

Symbols:
Qmin-Q25: from minimum read counts to 25% of the ROWS with most read counts
Q25-Q50:  from 25% read counts to 50% of the ROWS with most read counts
Q50-Q75:  from 50% read counts to 75% of the ROWS with most read counts
Q75-Qmax: from 75% read counts to maximum of the ROWS with most read counts ** top 25% suggested
Q50-Qmax: from 50% read counts to maximun of the ROWS with most read counts  ** top 50% suggested
Q25-Qmax: from 25% read counts to maxinum of the ROWS with most read counts  ** top 75% suggested
FR is the value 0/1 for either using filtering in ROWS or not.

FS is the value 0/1/2/3 for the (( filter_sStart )) input
FE is the value 0/1/2/3 for the (( filter_sEnd ))   input

Example 0: Qmin-Q25, figure: wochenende_heatmap_"number-label"_FR_FS_FE....0.png , values:
Example 1: Q25-Q50 , figure: wochenende_heatmap_"number-label"_FR_FS_FE....1.png , values:
Example 2: Q50-Q75 , figure: wochenende_heatmap_"number-label"_FR_FS_FE....2.png , values:
Example 3: Q75-Qmax, figure: wochenende_heatmap_"number-label"_FR_FS_FE....3.png , values:
Example 4: Q50-Qmax, figure: wochenende_heatmap_"number-label"_FR_FS_FE....4.png , values:
Example 5: Q25-Qmax, figure: wochenende_heatmap_"number-label"_FR_FS_FE....5.png , values:

argument: --filter_sTF      *TF = True or False *1 for True and 0 for False.
Should be 1 (one), If you are NOT willing to use all the files-samples in order to create heatmap, or those with only the files-samples with the most read counts.
Or should be 0 (zero)  If you are willing to use all the files-samples in order to create heatmap, or those with only the files-samples with the most read counts.

argument: --filter_sStart
Defines where the filtering should start.

argument: --filter_sEnd
Defines where the filtering should end.

argument: --suffix_label
Is a manual suffix to use into our .png file names just to avoid overwriting.

argument: --genus
0 for full name labeling of bacterial strain on each heatmap, 1 for Genus only per bacteria and 2 for Genus and Species per bacteria.

argument: --corr      *TF = True or False *1 for True and 0 for False.
If you wish to create simple heamap you should use 0 on this argument,
Or if you wish to create correlation heatmap you should use 1 on this argument.

argument: --where
In which directory does the user has his files -- where should python script look for the files.
________________________________________________________________________________

""",formatter_class=RawTextHelpFormatter)
##parser.add_argument('filename1',default='check_string_for_empty',
##                    help='After the python file .py you should write the name of the file of tsv format to read. This file should be in the same folder with the python.py file')
parser.add_argument('--batch_files', default="empty",
                     help='Please create a file that has in every line the filenames that you want to represent into the heatmap.This file should be an input to this argument.')
parser.add_argument('--filter_sTF', default=0,
                     help='Filter samples-patients yes:1 no:0 .Default: 0:No filtering')
parser.add_argument('--filter_sStart', default=0,
                     help='Filter samples-patients and where to start filtering from.Default:0 for Q0=min')
parser.add_argument('--filter_sEnd', default=4,
                     help='Filter samples-patients and whete to end filtering from.Default:4 for Q4=max')
parser.add_argument('--suffix_label', default=0,
                     help='Add any suffix label, used into the .png file names output. Default: 0')
parser.add_argument('--genus', default=0,
                     help='Labeling on heatmaps only with their Genus if input=1, or Genus_Species if input=2 or keep whole name if input=0. Default : 0')
parser.add_argument('--corr', default=0,
                     help='Make correlation heatmaps if input=1 or make simple heatmap sample-patients vs bacteria. Default : 0')
parser.add_argument('--where', default=os.getcwd(),
                     help="In which directory should we look for files.")

args = parser.parse_args()
#_____________________________________________________
##filename=args.filename1
fileloop=args.batch_files
filtering_samples=int(args.filter_sTF)
filtering_samplesEND=args.filter_sEnd
filtering_samplesSTART=args.filter_sStart
newfilenamesuffix=args.suffix_label
T_O_F_Genus=int(args.genus)
Corr_heat=int(args.corr)
wheredir=args.where


newfilenamesuffixor=newfilenamesuffix
#this is the suffix of the png creating heatmaps
#It includes the user given namefile, start quartile and then end quartile of samples-patients,
#1 for genus,0 for not genus labeling and 1 for correlation-heatmap, 0 for simple heatmaps
#and finally adds one final suffix of: 0-1-2-3 for 4 quartiles of totalrow counts per bacteria
newfilenamesuffix=str(str(newfilenamesuffix)+str(filtering_samplesSTART)
                    +str(filtering_samplesEND)+str(T_O_F_Genus)+str(Corr_heat))
filtering_samplesEND=int(filtering_samplesEND)+3     #describe() function input from 3=min to 7=max
filtering_samplesSTART=int(filtering_samplesSTART)+3 #describe() function input from 3=min to 7=max


#WHERE AM I???
#DO I HAVE A TSV FILE? #####IS IT A TSV?????############

dirpath = os.getcwd()
print("####Current directory is : " + dirpath)
print("####Current file.txt with filenames to open is:" + fileloop)

fileloop2=fileloop.split(".")[0]

try:
    fileloop2=fileloop.split(".")[0]
    fileloopJUSTCHECK=fileloop.split(".")[1]
except:
    fileloop2=str(fileloop[:4]+fileloop[-4:])
#############MAKE NEW DIRECTORIES

#-------------OUTPUT DIRECTORIES--------------

#THE HIGHEST 25% of totalread counts
def _25_of_all_(resultfilter2,describe_first,describe_second):
    Q1=resultfilter2['RowSum'].describe()[describe_first]
    Q2=resultfilter2['RowSum'].describe()[describe_second]
    higher_thenQ1=resultfilter2.loc[resultfilter2["RowSum"] >= Q1]
    higher_thenQ1=higher_thenQ1.loc[higher_thenQ1["RowSum"] <= Q2]
    higher_thenQ1=higher_thenQ1.drop("RowSum",1)
    return higher_thenQ1


I=["What"]
if fileloop=="empty":
    print("Something is wrong with the filenames, try to check if all filenames exist into your file.txt and if the named files exist")

else:
    f = open("{}/{}".format(dirpath,fileloop), "r")
    F=f.readlines()
    left=F[0].rstrip()
    left2=F[0].rstrip().split(".")[0]
    left3=left2.split("_")[0]+"_"+left2.split("_")[1]+"_"+left2.split("_")[2]
    print(left3,"dataset")               #for naming

    dataset_dir="{}/{}".format(dirpath,left3)
    try:
        os.mkdir(dataset_dir)
    except OSError:
        print ("Creation of the directory %s failed" % dataset_dir)
        print( '###__Maybe this directory already exists__###__OR you are not allowed to create it__###')
    else:
        print ("Successfully created the directory %s " % dataset_dir)

    dataset_heatmap_dir="{}/heatmap".format(dataset_dir)
    try:
        os.mkdir(dataset_heatmap_dir)
    except OSError:
        print ("Creation of the directory %s failed" % dataset_heatmap_dir)
        print( '###__Maybe this directory already exists__###__OR you are not allowed to create it__###')
    else:
        print ("Successfully created the directory %s " % dataset_heatmap_dir)

    leftPANDA=pd.read_csv("{}/{}".format(wheredir,left),delimiter="\t",header=None,usecols=[0,2],names=['colA', '{}'.format(left2)])
#        print(leftPANDA.head())
    for i in F[1:]:
        right=i.rstrip()
        right2=i.rstrip().split(".")[0]                #for labeling and heatmap figures
        rightPANDA=pd.read_csv("{}/{}".format(wheredir,right),delimiter="\t",header=None,usecols=[0,2],names=['colA', '{}'.format(right2)])
#            print(rightPANDA.head())
        I.append(i.rstrip())
        result = pd.merge(leftPANDA, rightPANDA, how='outer', on=['colA'])
        leftPANDA=result
    I=I[1:]
#        print(result)

    export_csv = result.to_csv ('{}/result.csv'.format(dataset_dir), index = None, header=True)
#    print(result)
result = result.set_index('colA')
result.index.names = [None]


#    To use the masked plot
result.fillna(0, inplace=True)
#+1 to all

listofChrom=['1_1_1_2', '1_1_1_1', '1_1_1_3', '1_1_1_4', '1_1_1_5', '1_1_1_6', '1_1_1_7', '1_1_1_8', '1_1_1_12', '1_1_1_11', '1_1_1_10', '1_1_1_9', '1_1_1_13','1_1_1_14', '1_1_1_18', '1_1_1_X', '1_1_1_15', '1_1_1_16', '1_1_1_17','1_1_1_20', '1_1_1_19', '1_1_1_21', '1_1_1_22', '1_1_1_Y', '1_1_1_MT']
resultfilter = result.drop(listofChrom, 0)
resultfilteradd1=resultfilter.add(1)
resultfilter = np.log2(resultfilteradd1)    # figure #2 was log 2 and figure #1$
print(resultfilter.max(),"this is the log2() ")
expone=2**resultfilter.max()
print(expone, "this is the Count")

result_mask = np.ma.masked_where(resultfilter == 0, resultfilter)
#    To use the unmasked plot--- the original heatmap

#empty for TOY DATASET###############################################################################
#if "MCF03s1_S8_R1" in list(resultfilter.columns.values):
#    resultfilter=resultfilter.drop("MCF03s1_S8_R1",1)
#######################################################################################################
print(resultfilter)
f.close()

plt.rcParams['figure.figsize'] = (16, 16)
my_dpi=96
multipdpi=4
plt.imshow(result_mask, interpolation = 'none', vmin = 0)
plt.savefig('{}/wochenende_{}_unmasked.png'.format(dataset_heatmap_dir,newfilenamesuffixor),dpi=my_dpi*multipdpi)
plt.close()
print(1*10 ," % Done")
#Z = linkage(resultfilter)  # You might want to set `method` and `metric`
#groups = fcluster(Z,t=0.8, criterion='distance')  # You might want to set `criterion`
#print(groups)

#np.savetxt("{}/groups2.txt".format(dataset_heatmap_dir),groups)

#groups = fcluster(Z,t=1, criterion='inconsistent')  # You might want to set `criterion`
#groups 1 = inconsistent
#groups 2 = distance
#list(resultfilter.rows.values)
#for future making new matrix col_A distance or cluster, col_B list of rows
col_list=list(resultfilter)
print(col_list)
#np.savetxt("{}/groups3.txt".format(dataset_heatmap_dir),groups)

##################ADD FILTER ON ROWS..... ->> so columns will be 50% not used.

if filtering_samples==1:
    resultfilter=resultfilter.T
    resultfilter["RowSum"]=resultfilter.sum(axis=1).astype(float)
    resultfilter.RowSum=resultfilter.RowSum.astype(float)
    resultfilter =_25_of_all_(resultfilter,filtering_samplesSTART,filtering_samplesEND)
    resultfilter=resultfilter.T

else:
   pass
####################################################
resultfilter["RowSum"]=resultfilter.sum(axis=1).astype(float)
resultfilter.RowSum=resultfilter.RowSum.astype(float)
##################
#GENUS????
if T_O_F_Genus==1 or T_O_F_Genus==2:
    resultfilter=resultfilter.reset_index()
    split_data = resultfilter["index"].str.split("_", n=5)
    data = split_data.to_list()
    names = ["NU","GENEBANKID","NU2","Genus","Species","What"]
    new_df = pd.DataFrame(data, columns=names)
    if T_O_F_Genus==2:
        new_df['Genus2'] = new_df[new_df.columns[3:5]].apply(
        lambda x: '_'.join(x.dropna().astype(str)),axis=1)
        resultfilter["Genus"]=new_df["Genus2"]
    else:
        resultfilter["Genus"]=new_df["Genus"]

    print(new_df)

    resultfilter.drop(columns =["index"], inplace = True)
    resultfilter=resultfilter.set_index("Genus")
else:
    pass

#####################################################

#####################
print(resultfilter,"AND")


import random
n1=len(set(list(resultfilter.index.values)))
def _colors_(n):
    ret = []
    N=n**(1./3.)
    #print(N, "cubic root for colors")
    r=0#1/N
    b=0#1/N
    g=0#1/N
    step = 1/(int(N)+1)######################################
    for i in range(int(N)+1):
        r=r+step
        for jj in range(int(N)+1):
            g=g+step
            for kk in range(int(N)+1):
                b=b+step

                ret.append((r,g,b))
            b=0
        g=0
    random.shuffle(ret)
    return ret
color1=_colors_(n1)
#print(color1)
#print(len(color1),"len of colors")

#print(color1)
#vmin = 4 :means 2 ^4 counts are visualized as not existing into our heatmap...
#Make a function-def
def _heatmap_(input_heatmap,name_heatmap,figure_size,xtick,ytick,cmap_color,multipdpi,my_dpi,T_O_F_Genus,Corr_heat):
    plt.rcParams['figure.figsize'] = (figure_size,figure_size)
    plt.rc('xtick', labelsize=xtick)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytick)    # fontsize of the tick labels
    fig,ax= plt.subplots()

    #input_heatmap2=input_heatmap.drop(columns=["Genus"])
    if Corr_heat==1:
        input_heatmap=input_heatmap.add(0.0001)
        input_heatmap=input_heatmap.T.corr()
        centerofcolor=int(0)
    else:
        centerofcolor=None

    if T_O_F_Genus==1:
        my_palette = dict(zip(input_heatmap.index.unique(), color1))
        row_colors =input_heatmap.index.map(my_palette)
        fig= sns.clustermap(input_heatmap,xticklabels=True, yticklabels=True,cmap=cmap_color
            ,row_colors=row_colors,center=centerofcolor)#,square=True)#,vmin=0) #metric="correlation", method="single"
    else:
        fig= sns.clustermap(input_heatmap,xticklabels=True, yticklabels=True,
                cmap=cmap_color,center=centerofcolor)#,square=True)#,vmin=0) #metric="correlation", method="single"

    plt.savefig('{}/wochenende_heatmap{}.png'.format(dataset_heatmap_dir,name_heatmap),dpi=my_dpi*multipdpi*2)
    plt.close("all")
    return


def _dendrodiagram_ (input_heatmap,name_heatmap,figure_size,xtick,ytick,multipdpi,my_dpi,T_O_F_Genus):
    plt.rcParams['figure.figsize'] = (figure_size,figure_size)

    Z = hierarchy.linkage(input_heatmap, 'ward')
    # Set the colour of the cluster here:
    #hierarchy.set_link_color_palette(['#b30000','#996600', '#b30086'])
    # Make the dendrogram and give the colour above threshold
    hierarchy.dendrogram(Z, color_threshold=240, above_threshold_color='grey',leaf_font_size=1.)
    # Add horizontal line.
    #plt.axhline(y=240, c='grey', lw=1, linestyle='dashed')
    plt.savefig('{}/wochenende_dendrogram{}.png'.format(dataset_heatmap_dir,name_heatmap),dpi=my_dpi*multipdpi*2)
    plt.close("all")
    return


if Corr_heat==1:
#factorfontsize=4/(filtering_samplesEND-filtering_samplesSTART)
    ffs=1
else:
    ffs=4/(filtering_samplesEND-filtering_samplesSTART)


###########################################################################
print(2*10 ," % Done")
newfilenamesuffix1=str(newfilenamesuffix)+str(0)
Qmin_Q25=_25_of_all_(resultfilter,3,4) # 3=min # 4=25pc # 5=50pc # 6=75pc # 7=max
#print(Qmin_Q25,"Qmin_Q25")
_heatmap_(Qmin_Q25,newfilenamesuffix1,3200,1.2*ffs,1.2,"coolwarm",4,96,T_O_F_Genus,Corr_heat)
print(3*10 ," % Done")
newfilenamesuffix1=str(newfilenamesuffix)+str(1)
Q25_Q50 =_25_of_all_(resultfilter,4,5) # 3=min # 4=25pc # 5=50pc # 6=75pc # 7=max
#print(Q25_Q50 , "Q25_Q50")
_heatmap_( Q25_Q50,newfilenamesuffix1,3200,1.2*ffs,1.2,"coolwarm",4,96,T_O_F_Genus,Corr_heat)
print(4*10 ," % Done")
newfilenamesuffix1=str(newfilenamesuffix)+str(2)
Q50_Q75 =_25_of_all_(resultfilter,5,6) # 3=min # 4=25pc # 5=50pc # 6=75pc # 7=max
#print(Q50_Q75 , "Q50_Q75")
_heatmap_( Q50_Q75,newfilenamesuffix1,3200,1.2*ffs,1.2,"coolwarm",4,96,T_O_F_Genus,Corr_heat)
print(5*10 ," % Done")
newfilenamesuffix1=str(newfilenamesuffix)+str(3)
Q75_Qmax=_25_of_all_(resultfilter,6,7) # 3=min # 4=25pc # 5=50pc # 6=75pc # 7=max
#print(Q75_Qmax,"Q75_Qmax")
_heatmap_(Q75_Qmax,newfilenamesuffix1,3200,1.2*ffs,1.2,"coolwarm",4,96,T_O_F_Genus,Corr_heat)
print(6*10 ," % Done")
newfilenamesuffix1=str(newfilenamesuffix)+str(4)
Q50_Qmax=_25_of_all_(resultfilter,5,7) # 3=min # 4=25pc # 5=50pc # 6=75pc # 7=max
#print(Q50_Qmax,"Q50_Qmax")
_heatmap_(Q50_Qmax,newfilenamesuffix1,200,1.2*ffs,1,"coolwarm",8,96,T_O_F_Genus,Corr_heat)
print(7*10 ," % Done")
newfilenamesuffix1=str(newfilenamesuffix)+str(5)
Q25_Qmax=_25_of_all_(resultfilter,4,7) # 3=min # 4=25pc # 5=50pc # 6=75pc # 7=max
#print(Q25_Qmax,"Q25_Qmax")
_heatmap_(Q25_Qmax,newfilenamesuffix1,200,1.2*ffs,1,"coolwarm",8,96,T_O_F_Genus,Corr_heat)
print(8*10 ," % Done")


####################################################################
#Still under construction:
######################################################################
newfilenamesuffix1=str(newfilenamesuffix)+str(6)
Z = hierarchy.linkage(Q75_Qmax, 'single')
plt.figure()
plt.rcParams['figure.figsize'] = (80,80)
#Now plot in given axes, improve the color scheme and use both vertical and horizontal orientations:
hierarchy.set_link_color_palette(['m', 'c', 'y', 'k','r','g','b','m', 'c', 'y', 'k','r','g','b'])
hierarchy.dendrogram(Z, above_threshold_color='y',
                           orientation='right',labels=Q75_Qmax.index.values,leaf_font_size=2)#,leaf_rotation=45)

hierarchy.set_link_color_palette(None)  # reset to default after use
my_dpi=96
multipdpi=1
plt.savefig('{}/wochenende_new_dendrogram_v3_{}.png'.format(dataset_heatmap_dir,newfilenamesuffix1),dpi=my_dpi*multipdpi*1)
plt.close()
print(9*10 ," % Done")

#Version 2:

# T. matrix for samples cluster-hierarchy-dendrogram
newfilenamesuffix1=str(newfilenamesuffix)+str(7)
Z = hierarchy.linkage(Q75_Qmax.T, 'single')
plt.figure()
plt.rcParams['figure.figsize'] = (200,200)
#Now plot in given axes, improve the color scheme and use both vertical and horizontal orientations:
hierarchy.set_link_color_palette(['m', 'c', 'y', 'k','r','g','b','m','#bcbddc' ,'c', 'y', 'k','r','g','b'])
hierarchy.dendrogram(Z, above_threshold_color='#bcbddc',
                           orientation='right',labels=Q75_Qmax.T.index.values)#,leaf_rotation=90)#,leaf_font_size=1.2)
hierarchy.set_link_color_palette(None)  # reset to default after use
my_dpi=96
multipdpi=1
plt.savefig('{}/wochenende_new_dendrogram_v3_{}.png'.format(dataset_heatmap_dir,newfilenamesuffix1),dpi=my_dpi*multipdpi*2)
plt.close()

#cmap_color GOOD: #BuPu # PiYG
print(10*10,"% DONE")
