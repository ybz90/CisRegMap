import numpy as np
import os
import subprocess
import sys
import time
from multiprocessing import Pool, Queue, Process, Manager, cpu_count
from scipy.stats import spearmanr as sp, ss



def pipe(cmd_list):
    return ' | '.join(cmd_list)

def run_bash(command):
    print command
    subprocess.call(command, shell=True)
    return None

def pretty(in_json):
    return json.dumps(in_json, indent=4, separators=(',', ': '))

def fast_pearson(curr_signal):

    # Calculate Pearson correlation pairwise over the rows (reg elements) of
    # the imported signal arrays
    start_time = time.time()

    # Define your input array of vectors
    #curr_signal = np.random.random((100, 1000))

    print curr_signal.shape
    rows = curr_signal.shape[0]

    # Init array to store pairwise correlation array
    pearson_array = []

    # Calculate Pearson row by row
    means = curr_signal.mean(axis=1)[(slice(None, None, None), None)]
    curr_signal_m = curr_signal - means
    curr_signal_ss = np.sqrt(ss(curr_signal_m, axis=1))

    # Write rows to pearson.txt; also export each row independently as its own
    # file, with row index # as filename
    for j in xrange(rows):
        # Calculate current row of Pearson correlations, p_row, to be written
        temp = np.dot(curr_signal_m, curr_signal_m[j].T)
        p_row = temp / (curr_signal_ss * curr_signal_ss[j])
        # print j, len(p_row)
        #p_row = np.around(p_row, decimals=10)
        pearson_array.append(p_row)

    pearson_array = np.array(pearson_array)
    print pearson_array.shape

    print ("--- %s seconds ---" % (time.time() - start_time))

    return pearson_array



###########
out_dir = '/Users/yuanz/Desktop/asdf'
cis_coords = '/Volumes/Seabiscuit/Work/Projects/ENCODE/ATAC-summits/all.tissues/summits/summits-merged.top.20.named.bed'
# NOTE: input coordinates must have name column; if len row == 3,
# append a name to each row.


'''Get the names for each element in coordinates list;
Save this dict to be used to retrieve row index later.'''

cis_names = {}
with open(cis_coords) as cis_coords_bed:
    for i, element in enumerate(cis_coords_bed):
        cis_names[element.rstrip().split('\t')[-1]] = i
# print np.array(cis_names)


'''Intersect coordinates of cis-regulatory elements with
interacting region boundaries to constrain interaction space.'''

# If using TAD domains...
boundary_coords = '/Users/yuanz/Desktop/asdf/total.HindIII.combined.domain'

# Otherwise, use arbitrary bins
bin_size = 500000
# From genome/chr_sz, generate bin coordinates, intersect as above

intersect_cmd = 'bedtools intersect -wo ' + \
    '-a %s ' % cis_coords + \
    '-b %s ' % boundary_coords + \
    '> %s/TAD_intersect.bed' % out_dir
run_bash(intersect_cmd)


'''Parse the intersect bed output to build dict of bin/TAD's and
the corresponding coordinates that fall within it.'''

TAD_loci_dict = {}

with open('%s/TAD_intersect.bed' % out_dir) as TAD_intersect:
    for i, line in enumerate(TAD_intersect):
        line = line.rstrip().split('\t')
        if '%s:%s-%s' % (line[-4], line[-3], line[-2]) not in TAD_loci_dict:
            TAD_loci_dict['%s:%s-%s' % (line[-4], line[-3], line[-2])] = []
        TAD_loci_dict['%s:%s-%s' %
                      (line[-4], line[-3], line[-2])].append(line[3])

print len(TAD_loci_dict)

# # debug -- see how many ele in each TAD; do name line up w/ idx
# z = [len(TAD_loci_dict[meep]) for meep in TAD_loci_dict]
# print np.array(z), np.max(z), np.min(z)

# zz = z.index(np.max(z))
# print zz

# asdf = [TAD_loci_dict[rawr] for rawr in TAD_loci_dict]
# flattened = [val for sublist in asdf for val in sublist]
# print np.array(flattened)
# print np.array([cis_names[asdf] for asdf in flattened])


'''Get the idx of "enhancer" and "promoter" elements relative
to the coordinates bed list.'''

# If using an input set of indices...
###

# Otherwise, define promoters as elements w/in 1kb of TSS
###

# Then, load the window bed output to get names of each se
enhancer_names = []
with open('%s/summits.named.TSS.distal.bed' % out_dir) as enhancer_bed:
    for j, line in enumerate(enhancer_bed):
        enhancer_names.append(line.rstrip().split('\t')[-1])
enhancer_idx = [cis_names[enhancer] for enhancer in enhancer_names]

promoter_names = []
with open('%s/summits.named.TSS.distal.bed' % out_dir) as promoter_bed:
    for j, line in enumerate(promoter_bed):
        promoter_names.append(line.rstrip().split('\t')[-1])
promoter_idx = [cis_names[promoter] for promoter in promoter_names]


'''Load signal array over which to run correlations.'''
# NOTE: In future versions, specify different arrays?
# NOTE2: Specify P-P, E-P, E-E, or all x all
# NOTE3: Choose Pearson or Spearman (or other?)

# Here, load signal array for ATAC-seq; get distal elements
# and proximal elements; correlate elements within same bin

ATAC_array = np.genfromtxt('%s/test_array.txt' % out_dir)
print ATAC_array.shape


'''For each bin/TAD, get the slice of input array(s) based on
corresponding row idx for constituent elements. Correlate the
elements of interest within the bin.'''

fast_pearson(np.random.random((1000, 66)))

TAD_corr = {}

for TAD in TAD_loci_dict:

    TAD_rows = TAD_loci_dict[TAD]
    TAD_rows_idx = [cis_names[ele_name] for ele_name in TAD_rows]

    TAD_array = np.array([ATAC_array[meep] for meep in TAD_rows_idx])

    TAD_correlation = fast_pearson(TAD_array)
    TAD_corr[TAD] = TAD_correlation


### NOTE: PLOT AND EXPORT SUMMARY DEBUG STATS/PLOTS HERE.


'''For promoter-enhancer predictions only!'''





