import matplotlib
matplotlib.use('Agg')

import glob
import os

import numpy as np
import seaborn as sns
import pandas as pd

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('mkrrgenome2bed.py -i [inputfile] -o [outputfile]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('mkrrgenome2bed.py -i [inputfile] -o [outputfile]')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print('Input file is:   ', inputfile)
    print('Output file is:  ', outputfile)

    df = pd.read_table(inputfile, header = None)

    def extract_start(x):
        x = x.split("..")[0]
        return x

    def extract_end(x):
        x = x.split(" ")[0]
        x = x.split("..")[1]
        return x

    def extract_end(x):
        x = x.split(" ")[0]
        x = x.split("..")[1]
        return x

    df[2] = df[1]
    df[1] = df[1].apply(lambda x: extract_start(x))
    df[2] = df[2].apply(lambda x: extract_end(x))

    df.to_csv(outputfile, header=None, index=None, sep='\t')

if __name__ == "__main__":
   main(sys.argv[1:])
