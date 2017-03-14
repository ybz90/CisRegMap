import json
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

# Calculate Pearson correlation pairwise over the rows (reg elements) of
# the imported signal arrays
def fast_pearson(curr_signal):

    #start_time = time.time()

    # Define your input array of vectors
    #curr_signal = np.random.random((100, 1000))

    #print curr_signal.shape
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
    #print pearson_array.shape

    #print ("--- %s seconds ---" % (time.time() - start_time))

    return pearson_array


###########

out_dir = '/Users/yuanz/Desktop/cisregmap'

cis_coords = '%s/summits-merged.top.20.named.bed' % out_dir

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
boundary_coords = '%s/mESC/HindIII_combined/total.HindIII.combined.domain' % out_dir

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

cis_types = {}

enhancer_names = []
with open('%s/summits-merged.top.20.named.distal.bed' % out_dir) as enhancer_bed:
    for j, line in enumerate(enhancer_bed):
        line = line.rstrip().split('\t')
        enhancer_names.append(line[-1])
        cis_types[line[-1]] = 'distal'

promoter_names = []
with open('%s/summits-merged.top.20.named.proximal.bed' % out_dir) as promoter_bed:
    for j, line in enumerate(promoter_bed):
        line = line.rstrip().split('\t')
        promoter_names.append(line[-1])
        cis_types[line[-1]] = 'proximal'

# print pretty(cis_types)


'''
Load signal array over which to run correlations.
Calculate pairwise correlation for EVERY locus within each TAD.
'''

# Here, load signal array for ATAC-seq; get distal elements
# and proximal elements; correlate elements within same bin

ATAC_array = np.genfromtxt('%s/test_array.txt' % out_dir)
print ATAC_array.shape


'''For each bin/TAD, get the slice of input array(s) based on
corresponding row idx for constituent elements. Correlate the
elements of interest within the bin.

Then, parse the pairwise correlation matrix for each interaction.
If the correlation is > threshold, add A->B and B->A to the
global interaction dict.'''

fast_pearson(np.random.random((1000, 66)))

r_cutoff = 0.7
interaction_pairs = {}

start_time = time.time()

# for TAD in ['chr1:3000000-4360000']:
for TAD in TAD_loci_dict:

    TAD_rows = TAD_loci_dict[TAD]

    if len(TAD_rows) > 1:
        TAD_rows_idx = [cis_names[ele_name] for ele_name in TAD_rows]

        TAD_array = np.array([ATAC_array[meep] for meep in TAD_rows_idx])

        # Calculate either Pearson or Spearman correlation
        # (or some alternative distance matrix)
        TAD_correlation = fast_pearson(TAD_array)

        # For each cell, the idx of curr TAD_array
        # corersponds to name in TAD_rows
        for i in range(len(TAD_rows)):
            for j in range(i+1, len(TAD_rows)):
                if TAD_correlation[i][j] > r_cutoff:  # add to interaction pairs dict
                    if TAD_rows[i] not in interaction_pairs:
                        interaction_pairs[TAD_rows[i]] = []
                    if TAD_rows[j] not in interaction_pairs:
                        interaction_pairs[TAD_rows[j]] = []
                    interaction_pairs[TAD_rows[i]].append(TAD_rows[j])
                    interaction_pairs[TAD_rows[j]].append(TAD_rows[i])

# print pretty(interaction_pairs)
print ('--- %s seconds ---' % (time.time() - start_time))


# NOTE: PLOT AND EXPORT SUMMARY DEBUG STATS/PLOTS HERE.

'''
For each interaction type, get the interactions from the global
interaction_pairs dict:
* E-P
* P-E
* P-P
* E-E
'''

EP = {}
PE = {}
PP = {}
EE = {}

for node in interaction_pairs:
    # Get the element type for each cis-element with predicted pairings
    node_type = cis_types[node]
    # Get the element type for each of its interacting targets
    node_targets = interaction_pairs[node]
    # Add the node, and its target interactions, to the appropriate dicts
    if node_type == 'distal':
        for target in node_targets:
            if cis_types[target] == 'distal':
                if node not in EE:
                    EE[node] = []
                EE[node].append(target)
            else:
                if node not in EP:
                    EP[node] = []
                EP[node].append(target)
    else:
        for target in node_targets:
            if cis_types[target] == 'distal':
                if node not in PE:
                    PE[node] = []
                PE[node].append(target)
            else:
                if node not in PP:
                    PP[node] = []
                PP[node].append(target)

#
print 'total interaction pairs: ', np.sum([len(interaction_pairs[node]) for node in interaction_pairs])
print 'num E-P: ', len(EP), len(EP) / float(len(enhancer_names)), np.sum([len(EP[node]) for node in EP])
print 'num P-E: ', len(PE), len(PE) / float(len(promoter_names)), np.sum([len(PE[node]) for node in PE])
print 'num P-P: ', len(PP), len(PP) / float(len(promoter_names)), np.sum([len(PP[node]) for node in PP])
print 'num E-E: ', len(EE), len(EE) / float(len(enhancer_names)), np.sum([len(EE[node]) for node in EE])


# NOTE: Export interaction allxall dict; desired interactions; plot a contact matrix like HiC.

