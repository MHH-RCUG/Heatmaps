# Heatmaps
Heatmap script for Wochenende output files containing abundance (reads per chromosome).

Input can be wochenende pipeline files (bam.txt files or filtered versions of these!).

## Installation

```
python3 -m pip install pandas
python3 -m pip install matplotlib
python3 -m pip install seaborn
```

## Simple usage

`
# Test run
python3 heat5.py --batch_files sample.txt --filter_sTF 1 --filter_sStart 0 --filter_sEnd 4 -- suffix_label A --genus 1 --corr 0
# Test run with example files
python3 heat5.py --batch_files eg1_sepsis_minion/sample.txt --filter_sTF 1 --filter_sStart 0 --filter_sEnd 4 -- suffix_label A --genus 1 --corr 0

`

## Tricks for heatmap (log2) etc.
To make it work i did those tricks into the python script:

1.  Kept 1 matrix with the "missing" counts, and made a figure out of that matrix. `**/wochenende_unmasked.png` 
2.  Replaced "missing" counts with zeros.
3.  Added 1 to all counts. (so that log2 of 1 would give me back 0 again)
4.  Calculated log2 to all counts (and replaced counts with their log2)
5.  Added 0.0001 to all 0 counts (otherwise heatmap would not work)
6.  (depending on the dataset and it's abundance of information)Filtered the samples having most read counts and the bacteria having most read counts too.

## Full usage:
```
Create heatmaps/correlation heatmaps and dendrograms generated from csv
files that are generated from wochenende pipeline. Those files should be filtered, sorted
 and be into a csv format as the above.

Example: MCF01_S7_R1.ndp.lc.trm.s.mq30.01mm.dup.calmd.bam.txt.filt.sort.csv

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Before running python script, make sure that you have a .txt file that incudes all the 
names-samples that you wish to combine, in order to create heatmaps.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
This .txt file should be used an argument --batch_files .txt (location)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   Each run of this script creates several heatmaps and dendrograms. 
In order to prevent overwriting, there are 2 different implemented label sets per .png files.

   1) A number-label that the user changes in every run that he/she decides.(input)
   2) A set of numbers that the original input for each argument user used or Default 
   values of the missing arguments, so that names of each heatmap can give information 
   for that heatmap.

    Each file has read counts into 1 column and those columns are merged into a matrix.
For instance, a matrix has rows that are bacterial strains and columns that are 
files-sample names.
The ammount of data that are used are divided into Quarters for rows and columns too!
        User can choose how many Quarters of the files-samples COLUMNS should be used 
for heatmap. 

Subsets can be set from argument (( filter_sStart )) and values according to limit are:
minimum:0
25%:1,
50%:2
75%:3
maximum:4

        User can choose any subset with those limits as starting poing
and as end point (( filter_sEnd ))
Of course it is suggested to have (( filter_sEnd )) > (( filter_sStart )) 
and (( filter_sEnd )) should be always set at value:4 in order to 
have the most informative files-samples

    Another filtering is done automatically and creates multiple 
heatmaps depending on the ROWS filtering per quartiles this time.
Different dimension on matrix for those Quarters then before, this 
time filtering is on ROWS and no input needed.

Symbols: 
RC is read counts.
Qmin-Q25: from minimum read counts to 25% of the ROWS with most read counts
Q25-Q50:  from 25% RC to 50% of the ROWS with most RC
Q50-Q75:  from 50% RC to 75% of the ROWS with most RC
Q75-Qmax: from 75% RC to maximum of the ROWS with most RC** top 25% suggested
Q50-Qmax: from 50% RC to maximun of the ROWS with most RC** top 50% suggested
Q25-Qmax: from 25% RC to maxinum of the ROWS with most RC** top 75% suggested
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
Should be 1 (one), If you are NOT willing to use all the files-samples 
in order to create heatmap, or those with only the files-samples with the most read counts.
Or should be 0 (zero)  If you are willing to use all the files-samples 
in order to create heatmap, or those with only the files-samples with the most read counts.

argument: --filter_sStart
Defines where the filtering should start.

argument: --filter_sEnd
Defines where the filtering should end.

argument: --suffix_label
Is a manual suffix to use into our .png file names just to avoid overwriting.

argument: --genus
0 for full name labeling of bacterial strain on each heatmap, 1 for Genus only 
per bacteria and 2 for Genus and Species per bacteria.

argument: --corr      *TF = True or False *1 for True and 0 for False.
If you wish to create simple heamap you should use 0 on this argument,
Or if you wish to create correlation heatmap you should use 1 on this argument.

argument: --where
In which directory does the user has his files -- where should python script 
look for the files.
```

## Examples and Disclaimers
```
Options:

Simple Heatmap (samples vs bacteria) or Correlation heatmap of bacteria vs bacteria ****
Filtered and kept 75% of top samples vs 25% of top bacteria or Filteted and kept 50% of top samples vs 25% of top bacteria.
Bacteria are displayed-labeled as Genus & Species or full name of genome strain.

wochenende_heatmap25024003.png      50% of top samples vs 25% of top bacteria , Full stain name,    simple Heatmap                             
wochenende_heatmap25024013.png      50% of top samples vs 25% of top bacteria , Full stain name,    Correlation Heatmap     
wochenende_heatmap25024203.png      50% of top samples vs 25% of top bacteria , Genus & Species, simple Heatmap      
wochenende_heatmap25024213.png      50% of top samples vs 25% of top bacteria , Genus & Species, Correlation Heatmap
wochenende_heatmap25014003.png      75% of top samples vs 25% of top bacteria , Full stain name,    simple Heatmap              
wochenende_heatmap25014013.png      75% of top samples vs 25% of top bacteria , Full stain name,    Correlation Heatmap
wochenende_heatmap25014203.png      75% of top samples vs 25% of top bacteria , Genus & Species, simple Heatmap    
wochenende_heatmap25014213.png      75% of top samples vs 25% of top bacteria , Genus & Species, Correlation Heatmap   




 +
wochenende_unmasked.png              ** 1 figure that can show the abundance of information into the original matrix.



Correlation matrix/heatmap:
****One issue about correlation is whether we accept that bacteria are coexisting with the
significant correlation that is displayed into the heatmaps, OR bacteria are "co-counted" 
with the same significant correlation because they share genes.
 
For example if 3-4 different bacteria of Streptococcus A-B-C-D species
appear into our sample then read counts of other Streptococcus E... species would be?? 
biased to have increased counts also, and not actually exist into our patient lungs..?



To conclude, this approach seems biased to the top read counts per bacteria and the 
top read counts per samples, since we are filtering out those that are not represented enough,
but it shows some clustering, clustering that was expected, and correlation between bacteria.

~Konstantinos~
```
