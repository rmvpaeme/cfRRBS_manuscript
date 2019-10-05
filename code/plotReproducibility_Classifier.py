## Import libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.switch_backend('Agg')
import glob
import os
import re
import sys, getopt
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess

### BEDtools > 2.27.1 is needed!!

# Folder to store intermediate files
tmp_folder = "./code/classifySamples/processed"

#sed -e "s/chr//" CpG_clusters.tsv > CpG_clusters_b37.tsv

clusters = pd.read_csv("./code/RRBS_450k_intersectClusters.tsv", sep="\t",usecols=[0,1,2], skiprows=[0], header=None, index_col=None)
clusters[3] = clusters.index
clusterFile = "RRBS_450k_intersectClusters"

clusters.to_csv(tmp_folder + "%s.bed" % clusterFile, header=None, index=None, sep='\t', mode = 'w')

# cfRRBS data
test_folder = "./code/classifySamples/testfiles"
test_files = glob.glob(os.path.join(test_folder, "*.cov"))
max = len(test_files)
count = 1

df_merged = pd.DataFrame()
for file in test_files:
        file_name = os.path.splitext(os.path.basename(file))[0]
        print("Importing %s" % file_name)
        print("File %i of %i (%2f)" % (count, max, (count/max)*100))
        df = pd.read_csv(file, sep="\t",usecols=[0,1,2,4,5], header=None, dtype =  {0: str, 1: np.str, 2: str} )
        df.to_csv(tmp_folder + "%s.bed" % file_name , header=None, index=None, sep='\t', mode = 'w')

        ## Get reads per clusters
        # Make a new file with the intersect of the clusters with the cov file
        outfile = open(tmp_folder + '%s_intersect.bed' % file_name, 'w')
        print("     Running bedtools intersect on %s.bed..." % file_name)
        arg = "bedtools intersect -wb -b %s/%s.bed -a %s%s.bed" % (tmp_folder, clusterFile, tmp_folder, file_name)
        arg = arg.split()
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        # The previous step shuffles the column order, so this step rearranges the column order
        df = pd.read_csv(tmp_folder + '%s_intersect.bed' % file_name, sep="\t", usecols=[6,7,8,3,4,5], header=None )
        df = df[[5,7,8,3,4,6]] # chr, start, stop, beta value, count methylated, count unmethylated, cluster number
        df.to_csv(tmp_folder + "%s_reordered.bed" % file_name , header=None, index=None, sep='\t', mode = 'w')

        # Group all the rows that are within a cluster, and get the sum of the methylated and unmethylated values. Also get the row with the CpG cluster number for future indexing.
        arg = "bedtools groupby -i %s%s_reordered.bed -g 1-3,6 -c 4,5 -o sum" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.bed' % file_name, 'w')
        print("     Running bedtools groupby on %s.bed..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_clustered.bed' % file_name, sep="\t", header=None, index_col = 2)

        df.index.name = None

        df = df[[0,1,3,4,5]] # chr, start, stop, beta value, no methylated, no unmethylated

        df.sort_values(by=[0,1,3], inplace=True)
        df.columns = ["Chr", "Start", "End", "Count Methyl", "Count Unmethyl"]
        df = df[pd.to_numeric(df["Chr"], errors = "coerce") < 23]
        df[file_name] = df["Count Methyl"] + df["Count Unmethyl"]
        df.sort_values(by=["Chr", "Start", "End"], inplace=True)
        df = df.drop(["Chr", "Start", "End", "Count Methyl", "Count Unmethyl"], axis = 1)
        if len(df_merged) == 0:
            df_merged = df
        else:
            df_merged = pd.merge(df, df_merged, how = "outer", left_index=True, right_index=True)
        count = count + 1

df_prct = pd.DataFrame()
df_prct["rco1"] = (df_merged.count(axis=1))
df_merged_mask30 = df_merged.mask(df_merged < 30)
df_prct["rco30"] = (df_merged_mask30.count(axis=1))
df_merged_mask100 = df_merged.mask(df_merged < 100)
df_prct["rco100"] = (df_merged_mask100.count(axis=1))
binsize = max + 1

countdf_rco1 = pd.DataFrame(df_prct["rco1"].value_counts(bins = range(0,binsize,1)))
countdf_rco1["bin"] = countdf_rco1.index.astype(str)
countdf_rco1["bin"] = countdf_rco1["bin"].apply(lambda x: x.split(',')[1])
countdf_rco1["bin"] =  countdf_rco1["bin"].apply(lambda x: x.split(']')[0]).astype(float)
countdf_rco1.sort_values(by=["bin"], inplace=True)
countdf_rco1.reset_index(inplace=True)
countdf_rco1.index = countdf_rco1["bin"]
countdf_rco1["rco1"] = countdf_rco1["rco1"].cumsum()
countdf_rco1["rco1"] = (countdf_rco1["rco1"]/countdf_rco1["rco1"].max())*100
countdf_rco1["bin"] = (countdf_rco1["bin"]/countdf_rco1["bin"].max())*100
countdf_rco1.drop(["index", "bin"], axis =1, inplace= True)

countdf_rco30 = pd.DataFrame(df_prct["rco30"].value_counts(bins = range(0,binsize,1)))
countdf_rco30["bin"] = countdf_rco30.index.astype(str)
countdf_rco30["bin"] = countdf_rco30["bin"].apply(lambda x: x.split(',')[1])
countdf_rco30["bin"] =  countdf_rco30["bin"].apply(lambda x: x.split(']')[0]).astype(float)
countdf_rco30.sort_values(by=["bin"], inplace=True)
countdf_rco30.reset_index(inplace=True)
countdf_rco30.index = countdf_rco30["bin"]
countdf_rco30["rco30"] = countdf_rco30["rco30"].cumsum()
countdf_rco30["rco30"] = (countdf_rco30["rco30"]/countdf_rco30["rco30"].max())*100
countdf_rco30["bin"] = (countdf_rco30["bin"]/countdf_rco30["bin"].max())*100
countdf_rco30.drop(["index", "bin"], axis =1, inplace= True)

countdf_rco100 = pd.DataFrame(df_prct["rco100"].value_counts(bins = range(0,binsize,1)))
countdf_rco100["bin"] = countdf_rco100.index.astype(str)
countdf_rco100["bin"] = countdf_rco100["bin"].apply(lambda x: x.split(',')[1])
countdf_rco100["bin"] =  countdf_rco100["bin"].apply(lambda x: x.split(']')[0]).astype(float)
countdf_rco100.sort_values(by=["bin"], inplace=True)
countdf_rco100.reset_index(inplace=True)
countdf_rco100.index = countdf_rco100["bin"]
countdf_rco100["rco100"] = countdf_rco100["rco100"].cumsum()
countdf_rco100["rco100"] = (countdf_rco100["rco100"]/countdf_rco100["rco100"].max())*100
countdf_rco100["bin"] = (countdf_rco100["bin"]/countdf_rco100["bin"].max())*100
countdf_rco100.drop(["index", "bin"], axis =1, inplace= True)

countdf = pd.merge(countdf_rco1, countdf_rco30, how = "inner", left_index=True, right_index=True)
countdf = pd.merge(countdf, countdf_rco100, how = "inner", left_index=True, right_index=True)
countdf["rco1"] =  abs(countdf["rco1"] - 100)
countdf["rco30"] =  abs(countdf["rco30"] - 100)
countdf["rco100"] =  abs(countdf["rco100"] - 100)
countdf["bin"] = countdf.index
countdf = countdf.melt('bin', var_name='rco',  value_name='count')
countdf = countdf.rename({"rco": "Read cutoff"}, axis='columns')
countdf.to_csv("reproducibility_matrix.csv", header=True, index = True, sep=',', mode = 'w')
