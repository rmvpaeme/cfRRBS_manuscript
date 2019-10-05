## Import libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import glob
import os
import re
import sys, getopt
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess
import argparse

# Parallelisation options
import multiprocessing
from multiprocessing import Process, Manager, Pool
cpuCount = (multiprocessing.cpu_count() - 2)

### BEDtools > 2.27.1 is needed!
### Tested with py >3.4
## More information is at https://github.com/shulik7/CancerLocator/tree/Add_instruction_and_data

## Read arguments from command line
parser = argparse.ArgumentParser('MakeTrainTest.py')
parser.add_argument('-t', '--train', default = 'disabled')
parser.add_argument('-v', '--viz', default = 'disabled')
args = parser.parse_args()
train = args.train
visualisation = args.viz

print("Generating the training set is:", train)

os.chdir("./")

# Name of output
testMethyName = "test_methy"
testDepthName = "test_depth"
trainFileName = "train"
testBetaName = "test_beta"

# Folder to store intermediate files
tmp_folder = "./classifySamples/processed/"

## The location of the cfRRBS files, after running the preprocessing pipeline
test_folder = "./classifySamples/testfiles/"
test_files = glob.glob(os.path.join(test_folder, "*.cov"))

# Infinium data for the reference dataset
NBL_infinium_folder = "./classifySamples/train/NBL"
NBL_infinium_files = glob.glob(os.path.join(NBL_infinium_folder, "*.txt"))
OS_infinium_folder = "./classifySamples/train/OS"
OS_infinium_files = glob.glob(os.path.join(OS_infinium_folder, "*.txt"))
CCSK_infinium_folder = "./classifySamples/train/CCSK"
CCSK_infinium_files = glob.glob(os.path.join(CCSK_infinium_folder, "*.txt"))
WBC_infinium_folder = "./classifySamples/train/WBC_child/"
WBC_infinium_files = glob.glob(os.path.join(WBC_infinium_folder, "*.txt"))
WT_infinium_folder = "./classifySamples/train/WT/"
WT_infinium_files = glob.glob(os.path.join(WT_infinium_folder, "*.txt"))
EWS_infinium_folder = "./classifySamples/train/EWS/450k"
EWS_infinium_files = glob.glob(os.path.join(EWS_infinium_folder, "*.txt"))
RMS_infinium_folder = "./classifySamples/train/RMS"
RMS_infinium_files = glob.glob(os.path.join(RMS_infinium_folder, "*.txt"))
MRT_infinium_folder = "./classifySamples/train/MRT"
MRT_infinium_files = glob.glob(os.path.join(MRT_infinium_folder, "*.txt"))

# WGBS data for the reference dataset.
NRML_WGBS_folder = "./classifySamples/train/NRML/"
NRML_WGBS_files = glob.glob(os.path.join(NRML_WGBS_folder, "*.cov"))

# This is a helper script, that if you have a m x n matrix (columns being samples and lines being CpGs), that splits the matrix into m x 1 and saves it back again. It is useful for preprocessing, but is not used here in the script.
# def splitMatrix(inputfile, outputfolder):
#     matrix = pd.read_csv(inputfile, sep=",", header=0, index_col=0)
#     for i in range(len(matrix.columns)):
#         name = list(matrix.columns.values)
#         name = name[i]
#         print("Writing %s.txt" % name)
#         matrix.to_csv(outputfolder + "%s.txt" % name , header=0, index=True,columns=[name], sep='\t', mode = 'w')

# The file containing the features (= the intersect between HM450K data and RRBS data, see GitHub README)
clusters = pd.read_csv("./classifySamples/resources/RRBS_450k_intersectClusters.tsv", sep="\t",usecols=[0,1,2], skiprows=[0], header=None, index_col=None)
clusters[3] = clusters.index
clusterFile = "RRBS_450k_intersectClusters"
clusters.to_csv(tmp_folder + "%s.txt" % clusterFile, header=None, index=None, sep='\t', mode = 'w')
clusters = clusters.drop([0,1,2,3], axis = 1) # Use empty index to later extract all the clusters from, so that every sample has the same number of clusters

# Load HumanMethylation450K reference file
array450k = pd.read_csv("./classifySamples/resources/HumanMethylation450_15017482_v1-2.csv", dtype={"CHR": str}, header = 7, usecols = (0,10,11,12), index_col="IlmnID")
array450k = array450k.dropna()
array450k[['MAPINFO', 'Genome_Build']] = array450k[['MAPINFO', 'Genome_Build']].astype(int)
array450k = array450k[array450k['Genome_Build'] == 37] # Extract locations with genome build GRCh37
array450k = array450k.drop(['Genome_Build'], axis = 1)
array450k[['CHR', 'MAPINFO']] = array450k[['CHR', 'MAPINFO']].astype(str)
array450k.index.name = None

# Load MethylationEPIC reference file
array850k = pd.read_csv("./classifySamples/resources/MethylationEPIC_v-1-0_B4.csv", dtype={"CHR": str}, header = 7, usecols = (0,10,11,12), index_col="IlmnID")
array850k = array850k.dropna()
array850k[['MAPINFO', 'Genome_Build']] = array850k[['MAPINFO', 'Genome_Build']].astype(int)
array850k = array850k[array850k['Genome_Build'] == 37] # Extract locations with genome build GRCh37
array850k = array850k.drop(['Genome_Build'], axis = 1)
array850k[['CHR', 'MAPINFO']] = array850k[['CHR', 'MAPINFO']].astype(str)
array850k.index.name = None


# Process test files for input in ruMeth_atlas.py
print("Generating %s and %s" % (testDepthName, testMethyName))
def import_test_files(x):
        # Goal: to obtain one file, containing all the test files, where the first column is all the samples and the rest of the column either the beta values, total depth or # methylated reads for that cluster.
        # 1. Read in the bismark coverage file and convert them to a sort-of bed file, so that it can be manipulated with bedtools.
        file = x
        file_name = os.path.splitext(os.path.basename(file))[0]
        df = pd.read_csv(file, sep="\t",usecols=[0,1,2,3,4,5], header=None)
        df[3] = df[3]/100    # From methylation percentage to methylation ratio
        df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

        # 2. Intersect between the cluster file and the bismark coverage file
        outfile = open(tmp_folder + '%s_intersect.txt' % file_name, 'w')
        print("     Running bedtools intersect on %s.txt..." % file_name)
        arg = "bedtools intersect -wb -b %s/%s.txt -a %s%s.txt" % (tmp_folder, clusterFile, tmp_folder, file_name)
        arg = arg.split()
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()
        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,7,8,3,4,5,9], header=None) # The previous step shuffles the column order, so this step rearranges the column order
        df = df[[6,7,8,3,4,5,9]] # chr, start, stop, beta value, count methylated, count unmethylated, cluster number
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

        # 3. Group all the rows that are within a cluster, and get the sum of the methylated and unmethylated values. In addition, get the row with the CpG cluster number for future indexing.
        arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,7 -c 5,6 -o sum" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
        print("     Running bedtools groupby on %s.txt..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        # 4. Remove all clusters that have less than 30 reads
        df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col = 3 )
        df.index.name = None    # Remove index.name for consistency
        df[6] = df[4]/(df[4] + df[5])   # Get beta value per cluster
        df = df[[0,1,2,6,4,5]] # Reorder the columns in chr, start, stop, beta value, no methylated, no unmethylated
        df.sort_values(by=[0,1,2], inplace=True) # Sort by chromosome
        df[7] = df[4] + df[5]   # Get total depth (=methylated + unmethylated count)
        df[[7,6,4,5]] = df[[7,6,4,5]].mask(df[7] < 30)  # Mark all clusters lower than 30 reads with NA
        print("The amount of clusters in %s remaining after removing NA: %s" % (file_name, len(df[7].replace(np.nan, 'NA').apply(pd.to_numeric, errors='coerce').dropna())))
        df = df.replace(np.nan, 'NA', regex=True) # Replace numpy NaN with NA string
        print("     Extracting %s and %s from file %s.txt..." % (testMethyName, testDepthName, file_name))

        # 5. Add the first file to the testMethy_list ( = number of methylated CpGs per cluster)
        testMethy_df = df
        testMethy_df.columns = [0,1,2,6, "%s" % file_name, 5,7]
        testMethy_df = testMethy_df.drop([0,1,2,6,5,7], axis = 1).astype(str)
        testMethy_df[file_name] = testMethy_df[file_name].apply(lambda x: x.split('.')[0]) # Make integer of float numbers, but as the dataframe is astype(str), we need to do this with a lambda function
        testMethy_list.append(testMethy_df)

        # Identical to testMethy
        testDepth_df = df
        testDepth_df.columns = [0,1,2,5,3,4,"%s" % file_name]
        testDepth_df = testDepth_df.drop([0,1,2,3,4,5], axis = 1).astype(str)
        testDepth_df[file_name] = testDepth_df[file_name].apply(lambda x: x.split('.')[0])
        testDepth_list.append(testDepth_df)

        # Make a new variable for visualisation that contains the beta values of the clusters
        testBeta_df = df
        testBeta_df.columns = [0,1,2, "%s" % file_name, 4, 5,7]
        testBeta_df = testBeta_df.drop([0,1,2,4,5,7], axis = 1).astype(str)
        testBeta_list.append(testBeta_df)

# Use the manager package so the lists are shared between the processes
with Manager() as manager:
    # Define empty lists
    testMethy_list = manager.list()
    testDepth_list = manager.list()
    testBeta_list  = manager.list()

    pool = Pool(cpuCount)  # Parallelisation function

    pool.map(import_test_files, test_files)    # Import the files in parallel

    print("Merging all %s in one file..." % testMethyName)  # Merge the testMethy_list in a pandas dataframe, merging the same indices
    testMethy = pd.concat(testMethy_list, axis = 1)
    testMethy = pd.merge(clusters, testMethy, how = "left", left_index=True, right_index=True)    # Merge the pandas df with the clusters, leaving NA values for clusters that were not covered.
    testMethy = testMethy.transpose().fillna('NA')   # Transpose and Fill NaN with NA string
    testMethy.to_csv("./classifySamples/output/%s" % testMethyName, header=None,sep='\t', mode = 'w')

    print("Merging all test_depth in one file...")
    testDepth = pd.concat(testDepth_list, axis = 1)
    testDepth = pd.merge(clusters, testDepth, how = "left", left_index=True, right_index=True)
    testDepth = testDepth.transpose().fillna('NA')
    testDepth.to_csv("./classifySamples/output/%s" % testDepthName, header=None,sep='\t', mode = 'w')

    testBeta = pd.concat(testBeta_list, axis = 1)
    testBeta = pd.merge(clusters, testBeta, how = "left", left_index=True, right_index=True)
    testBeta = testBeta.transpose().fillna('NA')
    print("Writing to disk...")
    testBeta.to_csv("./classifySamples/output/%s" % testBetaName, header=None,sep='\t', mode = 'w')

    testMethy_rmNA = testMethy.apply(pd.to_numeric, errors='coerce').dropna(axis=1)
    testDepth_rmNA = testDepth.apply(pd.to_numeric, errors='coerce').dropna(axis=1)
    print("The number of columns in the %s file after removing NA values is: %i" % (testMethyName,testMethy_rmNA.shape[1]))
    print("The number of columns in the %s file after removing NA values is: %i" % (testDepthName,testDepth_rmNA.shape[1]))

# Make the reference set, only needs to be done once, or every time new files are added to the reference dataset. Because the reference set is derived from multiple sources of publicly available data, every reference entity has it's own import (and no universal function for was written, because that would need additional preprocessing steps).
if train == 'enabled':
    print("Generating %s" % trainFileName)
    def getAvg(x):
        ## This function gets the average beta value in a cluster, or writes NA if more than half are not available.
        if isinstance(x, float):
            return x
        elif len(x) == 0:
            x = "NA"
        else:
            line_values = []
            countNA = 0
            line_values = x.split(',')
            line_values = [i.strip(' ') for i in line_values]
            line_values = list(filter(None, line_values))
            countTot = len(line_values)
            for value in line_values:
                if value == 'NA':
                    countNA = countNA +1
            if countTot == 0:
                x = "NA"
            elif countNA/countTot >= 0.5:
                x = "NA"
            else:
                line_values_rmNA = filter(lambda a: a != 'NA', line_values)
                line_values_rmNA_list = list(line_values_rmNA)
                calcmean = np.array(line_values_rmNA_list).astype(np.float)
                x = np.mean(calcmean)
            return x


    def generateTrain_Infinium(label, file_name):
        # The input for this function is an ordered infinium 450k file with the order chr - start - stop - beta value. If it doesn't have this structure, some preprocessing needs to be done.
        outfile = open(tmp_folder + "%s_intersect.txt" % file_name, 'w')
        print("     Running bedtools intersect on %s.txt..." % file_name)
        proc = subprocess.Popen(args=["bedtools", "intersect", "-b", tmp_folder + "%s.txt" % clusterFile, "-a", tmp_folder + "%s.txt" % file_name, "-wb"], stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,3,4,5,7], header=None )
        df = df[[4,5,6,3,7]]
        df[3] = df[3].replace(np.nan, 'NA', regex=True)
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

        # Group the total and methylated reads per cluster
        arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,5 -c 4 -o collapse" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
        print("     Running bedtools groupby on %s.txt..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col=3 )
        df[4] = df[4].astype(str)
        df = df.groupby([df.index,0,1,2])[4].apply(','.join) # BEDtools groupby doesnt always make the index unique for some reason, this fixes this.
        df = df.reset_index().set_index(3)
        df = df[~df.index.duplicated(keep='first')] # This shouldn't be necessary anymore, but keep it here as a double check
        df.index.name = None
        df.sort_values(by=[0,1,2], inplace=True)
        df[4] = df[4].apply(getAvg)
        df.columns = [0,1,2,label]
        df = df.drop([0,1,2], axis = 1)
        return df

    def generateTrain_NGS(inputfile, label, file_name):
        # The structure of this function is very similar to import_test_files()
        df = pd.read_csv(inputfile, sep="\t",usecols=[0,1,2,3,4,5], header=None)
        df[3] = df[3]/100
        df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

        outfile = open(tmp_folder + '%s_intersect.txt' % file_name, 'w')
        print("     Running bedtools intersect on %s.txt..." % file_name)
        arg = "bedtools intersect -wb -b %s%s.txt -a %s%s.txt" % (tmp_folder, clusterFile, tmp_folder, file_name)
        arg = arg.split()
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,7,8,3,4,5,9], header=None )
        df = df[[6,7,8,3,4,5,9]]
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name, header=None, index=None, sep='\t', mode = 'w')

        arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,7 -c 5,6 -o sum" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
        print("     Running bedtools groupby on %s.txt..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col = 3 )
        df.index.name = None

        df[6] = df[4]/(df[4] + df[5])  # Get beta value per cluster

        df = df[[0,1,2,6,4,5]]         # Reorder
        df.sort_values(by=[0,1,2], inplace=True)

        df[7] = df[4] + df[5]         # Get total depth

        df[[7,6,4,5]] = df[[7,6,4,5]].mask(df[7] < 30)         # Mark all clusters lower than 30 reads with NA
        df = df.replace(np.nan, 'NA', regex=True)
        df.columns = [0,1,2,label,4,5,7]
        df = df.drop([0,1,2,4,5,7], axis = 1)
        return df

    # Similar to import_test_files
    with Manager() as manager:
        #Define empty list
        trainFile_list = manager.list()

        # Read in WBC files.
        def import_WBC_train(x):
            # 1. First, some reordering is done so that the files can be manipulated with bedtools.
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Open the file
            df = pd.read_csv(file, sep="\t", header=None, index_col=0, names = ["Beta_Value"])
            ## Add the chromosomal position to the sample
            df = pd.merge(array450k, df, how = "inner", left_index=True, right_index=True)
            ## Add a stop and reorder the columns
            df["MAPINFO_Stop"] = df["MAPINFO"]
            df = df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_Value"]]
            df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            # 2. See function above, groups the cg-sites into clusters and merges it to a list with the prespecified label.
            df = generateTrain_Infinium(label = "wbc", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_WBC_train, WBC_infinium_files)

        ### Read in NBL files
        ## Specify files with clinical parameters (MYCN amplified or non-amplified, each line is a sample name)
        MNA = pd.read_table('./classifySamples/train/MNA.txt', sep='\n', header = None)
        MA = pd.read_table('./classifySamples/train/MA.txt', header = None)

        def import_NBL_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            name = file_name.split('.')
            name = str(name[4])
            name = "-".join(name.split("-", 3)[:3])
            # Extract position and beta value
            df = pd.read_csv(file, sep="\t",usecols=[0,1,2,3,4], header=1)
            df["Genomic_Coordinate_Stop"] = df["Genomic_Coordinate"]
            df = df[["Chromosome", "Genomic_Coordinate", "Genomic_Coordinate_Stop", "Beta_value"]]
            df.sort_values(by = ["Chromosome", "Genomic_Coordinate"], inplace=True)
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            if name in MNA.values:
                df = generateTrain_Infinium(label = "nbl-mna", file_name = file_name)
                trainFile_list.append(df)
            elif name in MA.values:
                df = generateTrain_Infinium(label = "nbl-ma", file_name = file_name)
                trainFile_list.append(df)
            else:
                df = generateTrain_Infinium(label = "nbl-nos", file_name = file_name)
                trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_NBL_train, NBL_infinium_files)

        # Similar to NBL, only the column headers are different.
        def import_OS_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Extract position and beta value
            df = pd.read_csv(file, sep="\t",usecols=[1,2,3,4,5], header=0)
            ## Add a stop and reorder the columns
            df["Position_Stop"] = df["Position"]
            df = df[["Chromosome", "Position", "Position_Stop", "Signal"]]
            df.sort_values(by = ["Chromosome", "Position"], inplace=True)
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            df = generateTrain_Infinium(label = "os", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_OS_train, OS_infinium_files)

        # # Similar to NBL, only the column headers are different.
        def import_CCSK_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Extract position and beta value
            df = pd.read_csv(file, sep="\t",usecols=[1,2,3,4,5], header=0, dtype={"Position": str})
            df
            ## Add a stop and reorder the columns
            df["Position_Stop"] = df["Position"]
            df = df[["Chromosome", "Position", "Position_Stop", "AVG_Beta"]]
            df.sort_values(by = ["Chromosome", "Position"], inplace=True)
            df = df.dropna()
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            df = generateTrain_Infinium(label = "ccsk", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_CCSK_train, CCSK_infinium_files)

        ### Similar to NBL, only the column headers are different.
        def import_WT_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Extract position and beta value
            df = pd.read_csv(file, sep="\t",usecols=[1,2,3,4,5], header=0, dtype={"Position": str, "Chromosome": str})
            ## Add a stop and reorder the columns
            df["Position_Stop"] = df["Position"]
            df = df[["Chromosome", "Position", "Position_Stop", "AVG_Beta"]]
            df.sort_values(by = ["Chromosome", "Position"], inplace=True)
            df = df.dropna()
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            df = generateTrain_Infinium(label = "wt", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_WT_train, WT_infinium_files)

        def import_EWS_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Open the file
            df = pd.read_csv(file, sep="\t", header=None, index_col=0, names = ["Beta_Value"])
            ## Add the chromosomal position to the sample
            df = pd.merge(array450k, df, how = "inner", left_index=True, right_index=True)
            ## Add a stop and reorder the columns
            df["MAPINFO_Stop"] = df["MAPINFO"]
            df = df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_Value"]]
            df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            df = generateTrain_Infinium(label = "ews", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_EWS_train, EWS_infinium_files)

        def import_MRT_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Open the file
            df = pd.read_csv(file, sep="\t", header=None, index_col=0, names = ["Beta_Value"])
            ## Add the chromosomal position to the sample
            df = pd.merge(array850k, df, how = "inner", left_index=True, right_index=True)
            ## Add a stop and reorder the columns
            df["MAPINFO_Stop"] = df["MAPINFO"]
            df = df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_Value"]]
            df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            df = generateTrain_Infinium(label = "mrt", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_MRT_train, MRT_infinium_files)

        ## Specify files with clinical parameters
        aRMS = pd.read_table('./classifySamples/train/aRMS.txt', sep='\n', header = None)
        eRMS = pd.read_table('./classifySamples/train/eRMS.txt', header = None)

        def import_RMS_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            ## Open the file
            df = pd.read_csv(file, sep="\t", header=None, index_col=0, names = ["Beta_Value"])
            ## Add the chromosomal position to the sample
            df = pd.merge(array450k, df, how = "inner", left_index=True, right_index=True)
            ## Add a stop and reorder the columns
            df["MAPINFO_Stop"] = df["MAPINFO"]
            df = df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_Value"]]
            df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
            df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

            if file_name in aRMS.values:
                df = generateTrain_Infinium(label = "arms", file_name = file_name)
                trainFile_list.append(df)
            elif file_name in eRMS.values:
                df = generateTrain_Infinium(label = "erms", file_name = file_name)
                trainFile_list.append(df)
            else:
                df = generateTrain_Infinium(label = "rms-nos", file_name = file_name)
                trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_RMS_train, RMS_infinium_files)

        def import_NRML_train(x):
            file = x
            file_name = os.path.splitext(os.path.basename(file))[0]
            df = generateTrain_NGS(inputfile = file, label = "normal", file_name = file_name)
            trainFile_list.append(df)

        pool = Pool(cpuCount)
        pool.map(import_NRML_train, NRML_WGBS_files)

        # Generate full matrix from list
        trainFile = pd.concat(trainFile_list, axis = 1)
        # Make sure that the file contains all the clusters
        trainFile = pd.merge(clusters, trainFile, how = "left", left_index=True, right_index=True)
        trainFile = trainFile.transpose().fillna('NA')
        trainFile.to_csv("./classifySamples/output/%s" % trainFileName, header=None, sep='\t', mode = 'w')
        trainFile_rmNA = trainFile.select_dtypes(include=['float64'])
        print("The number of columns in the %sfile is: %i" %  (trainFileName,trainFile.shape[1]))
        print("The number of columns in the %sfile after removing all NA values is: %i" %  (trainFileName,trainFile_rmNA.shape[1]))

    if visualisation == 'enabled':
        trainFile = trainFile.transpose()
        testBeta = testBeta.transpose()
        TotalMatrix = trainFile

        outputfolder = "./classifySamples/output/plots"
        TotalMatrix = TotalMatrix.apply(pd.to_numeric, errors='coerce').dropna()
        print("The number of remaining rows in the clustermap is %d (after removing rows containing NA values)" % (len(TotalMatrix)))
        TotalMatrix.to_csv('%s/TotalMatrix.csv' % outputfolder, sep=',', mode='w')

        # Get colors for each tumor in the plot
        TotalMatrix_labels = TotalMatrix.columns.unique()
        TotalMatrix_pal = sns.cubehelix_palette(TotalMatrix_labels.unique().size,
                                            light=.9, dark=.1, reverse=True,
                                            start=1, rot=-2)
        TotalMatrix_lut = dict(zip(map(str, TotalMatrix_labels.unique()), TotalMatrix_pal))
        TotalMatrix_colors = pd.Series(TotalMatrix_lut)


        print("Generating tSNE")
        TotalMatrix = TotalMatrix.transpose()
        ## Make a new column with the name of the indices
        TotalMatrix['index1'] = TotalMatrix.index
        ## Extract the tumor name from the indices
        TotalMatrix['tumor'] = TotalMatrix['index1']
        TotalMatrix.drop('index1', axis = 1, inplace = True)

        matplotlib.rcParams.update({'font.size': 18})
        ### tSNE plots
        X_tsne = TotalMatrix.drop("tumor", axis = 1)
        y_tsne = TotalMatrix['tumor']
        print("The number of CpGs in the tSNE-plot is: %i" % len(TotalMatrix.columns))
        from sklearn.manifold import TSNE
        tsne = TSNE(n_components=2, verbose=1, perplexity=30, n_iter=2000)
        tsne_results = tsne.fit_transform(X_tsne)
        tsneDf = pd.DataFrame(data = tsne_results
                     , columns = ['t-SNE 1', 't-SNE 2'])
        y_tsne = y_tsne.values
        final_tSNE_Df = pd.concat([tsneDf, pd.DataFrame(y_tsne, columns=["target"])], axis = 1)
        print("Generating t-SNE plots...")
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('t-SNE 1', fontsize = 18)
        ax.set_ylabel('t-SNE 2', fontsize = 18)
        ax.set_title('2 component t-SNE', fontsize = 20)
        targets = list(TotalMatrix['tumor'].unique())
        colors = TotalMatrix_pal
        for target, color in zip(targets,colors):
            indicesToKeep = final_tSNE_Df['target'] == target
            print('Number of %s on t-SNE plot: %d' % (target,len(final_tSNE_Df.loc[indicesToKeep])))
            ax.scatter(final_tSNE_Df.loc[indicesToKeep, 't-SNE 1']
                       , final_tSNE_Df.loc[indicesToKeep, 't-SNE 2']
                       , c = color #pd.Series({"wbc": "#43b7ba", "nbl-ma": "#ec672e", "nbl-mna": "#1b2944", "NBL1-cfRRBS":"#000000", "NBL2-cfRRBS":"#000000", "NBL1-WGBS":"#000000", "NBL2-WGBS":"#000000", "NBL1-SeqCapEpi":"#000000", "NBL2-SeqCapEpi":"#000000"})
                       , s = 30)
        ax.legend(targets, fontsize = 16, bbox_to_anchor=(1.04,1), loc="upper left", ncol=1, fancybox=True)
        ax.grid(False)
        fig.show()
        fig.savefig('%s/tSNEplot.png' % outputfolder, dpi = 300, bbox_inches="tight")
        fig.savefig('%s/tSNEplot.svg' % outputfolder, dpi = 300, bbox_inches="tight")
        print("Done with plotting tSNE.")


        import umap
        X_umap = TotalMatrix.drop("tumor", axis = 1)
        y_umap = TotalMatrix['tumor']

        print("The number of CpGs in the UMAP-plot is: %i" % len(TotalMatrix.columns))
        umap_results = umap.UMAP().fit_transform(X_umap)

        umapDf = pd.DataFrame(data = umap_results
                     , columns = ['UMAP 1', 'UMAP 2'])
        y_umap = y_umap.values
        final_umap_Df = pd.concat([umapDf, pd.DataFrame(y_umap, columns=["target"])], axis = 1)
        print("Generating t-SNE plots...")
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('UMAP 1', fontsize = 18)
        ax.set_ylabel('UMAP 1', fontsize = 18)
        ax.set_title('2-component UMAP', fontsize = 20)
        targets = list(TotalMatrix['tumor'].unique())
        colors = TotalMatrix_pal
        for target, color in zip(targets,colors):
            indicesToKeep = final_umap_Df['target'] == target
            print('Number of %s on UMAP plot: %d' % (target,len(final_umap_Df.loc[indicesToKeep])))
            ax.scatter(final_umap_Df.loc[indicesToKeep, 'UMAP 1']
                       , final_umap_Df.loc[indicesToKeep, 'UMAP 2']
                       , c = color
                       , s = 30
                       )
        #ax.legend(targets, fontsize = 16)
        ax.legend(targets, fontsize = 16, bbox_to_anchor=(1.04,1), loc="upper left",
              ncol=1, fancybox=True)
        ax.grid(False)
        fig.show()
        fig.savefig('%s/UMAPplot.svg' % outputfolder, dpi = 300, bbox_inches="tight")
        fig.savefig('%s/UMAPplot.png' % outputfolder, dpi = 300, bbox_inches="tight")
        print("Done with plotting UMAP.")

print("The number of clusters is: %i" % len(clusters))
print("Done")
