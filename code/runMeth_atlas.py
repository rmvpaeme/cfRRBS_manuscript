## Import libraries
import glob
import os
import re
import numpy as np
import seaborn as sns
import pandas as pd
import subprocess
import deconvolve_resids

## Set the working directory
os.chdir('/Users/rmvpaeme/Repos/cfRRBS_paper/code/')

# Input files (= outputfiles from MakeTrainTest.py)
inputFile_tumor_atlas = "./classifySamples/output/train"

inputFile_samples = "./classifySamples/output/test_beta"
# inputFile_samples = "./classifySamples/output/test_beta_manuscript" # to reproduce manuscript results

outputFile_atlas = "train_methatlas.csv"
outputFile_samples = "test_methatlas.csv"
outDir = "./classifySamples/output/classification"

classificationResults = "classificationResults_methAtlas.csv"

df_tumor = pd.read_csv(inputFile_tumor_atlas, sep="\t", header = None, index_col = None)
df_tumor = df_tumor.transpose()
df_tumor.columns = df_tumor.iloc[0]
df_tumor = df_tumor.reindex(df_tumor.index.drop(0))
df_tumor = df_tumor.astype('float64')
df_tumor = df_tumor.groupby(by=df_tumor.columns, axis=1).median()
df_tumor['IlmnID'] = df_tumor.index
df_tumor['IlmnID'] = 'cg' + df_tumor['IlmnID'].astype(str)
df_tumor = df_tumor.set_index("IlmnID")

df_tumor.to_csv("%s/%s" % (outDir,outputFile_atlas), header=True, index = True, sep=',', mode = 'w')
# df_tumor = pd.read_csv("./train_methatlas_plasma_manuscript.csv", sep=",", header = 0, index_col = 0) # uncomment this to reproduce the manuscript results

df_samples = pd.read_csv(inputFile_samples, sep="\t", header = None, index_col = None)
df_samples = df_samples.transpose()
df_samples.columns = df_samples.iloc[0]
df_samples = df_samples.reindex(df_samples.index.drop(0))
df_samples = df_samples.astype('float64')
df_samples['IlmnID'] = df_samples.index
df_samples['IlmnID'] = 'cg' + df_samples['IlmnID'].astype(str)
df_samples = df_samples.set_index("IlmnID")
df_samples = df_samples.reindex(sorted(df_samples.columns), axis=1)
df_samples.to_csv("%s/%s" % (outDir,outputFile_samples), header=True, index = True, sep=',', mode = 'w')

deconvolve_resids.Deconvolve(atlas_path="%s/%s" % (outDir,outputFile_atlas),samp_path="%s/%s" % (outDir,outputFile_samples),out_dir=outDir,plot=False,resid=True).run()

results = pd.read_csv("%s/%s" % (outDir,outputFile_samples.split('.')[0] + "_deconv_output.csv"), sep=",", header = 0, index_col = 0)
results = results.drop(["wbc", "normal"])
results = results.replace(0, np.nan)
results.sort_values(by=list(results.columns.values), inplace=True)

dict = {}
for column in results:
    dict[column] = [results[column].idxmax(), results[column].max()]
    results_tumor = pd.DataFrame.from_dict(dict, orient = "index")
results_tumor = results_tumor.rename(columns={0: "Classification", 1: "Tumor burden"})
results_tumor = results_tumor.replace(np.nan, 0)
results_tumor = results_tumor.sort_index()
results_tumor.to_csv("%s/%s" % (outDir,classificationResults), header=True, index = True, sep=',', mode = 'w')
