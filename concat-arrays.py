# Import modules

import csv
import fastcluster
import h5py
import itertools
import json
import matplotlib
matplotlib.use('pdf')
%matplotlib inline
import numpy as np
import os
import pandas as pd
import pylab
import scipy.cluster.hierarchy as sch
import seaborn as sns
import subprocess
import sys
import time
from IPython.display import clear_output
from matplotlib import pyplot as plt
from multiprocessing import Pool, Queue, Process, Manager, cpu_count
from natsort import natsorted
from numpy import random
from pandas import DataFrame as df
from scipy.stats import pearsonr as pc, spearmanr as sp, norm, rayleigh, gaussian_kde, rankdata, ss
from scipy.cluster.vq import *

print sys.getrecursionlimit()
sys.setrecursionlimit(20000)

#
def pipe(cmd_list):
    return ' | '.join(cmd_list)

#
def run_bash(command):
    print command
    subprocess.call(command, shell=True)
    return None
#
def pretty(in_json):
    return json.dumps(in_json, indent=4, separators=(',', ': '))


'''
Initialize global variables and parameters.
'''

# Project root directory
proj_root = '/Volumes/Seabiscuit/Work/ENCODE/'

# Define current analysis directory; get tissues in directory
out_root = proj_root + 'ATAC-summits'
tissues = [tissue for tissue in os.listdir(
    out_root) if '.DS_Store' and '.' not in tissue]
print tissues

# Define histone marks (or other datasets e.g. RNAseq)
#marks = ['ATAC','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','H3K27me3','H3K36me3']
marks = ['ATAC']


'''Scan each tissue /peaks folder for stages to parse, cat the pooled.truerep
narrowPeak predictions, and then sort/merge/name the peaks to produce the
uniform, merged peak set.'''

tissue_stages = {}
for tissue in tissues:
    stages = [stage for stage in os.listdir('%s/%s/peaks' % (out_root, tissue))
              if stage != '.DS_Store'
              if '.bed' not in stage]
    tissue_stages[tissue] = stages

print pretty(tissue_stages)



##
# SPECIFIC ORDER ECTO, MESO, ENDO -- NOT JUST ALPHA
tissues2 = ['forebrain','midbrain','hindbrain','neural tube',
            'limb','embryonic facial prominence','heart','liver',
            'intestine','kidney','lung','stomach']


# load the entire catalog array; for each tissue and cat it all together

all_stacked = np.zeros((336484,0))

for tissue in tissue_stages:

    tissue_stacked = np.genfromtxt('%s/%s/arrays/ATAC/mean.quantile.top.50.txt' %(out_root, tissue))
    print tissue, tissue_stacked.shape

    all_stacked = np.hstack([all_stacked, tissue_stacked])

print all_stacked.shape


np.savetxt('/Users/yuanz/Desktop/test_array.txt',
       all_stacked, fmt='%1.6f', delimiter='\t')

