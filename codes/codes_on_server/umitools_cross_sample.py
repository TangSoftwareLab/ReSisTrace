############################################################
# Using umi_tools to merge lineage label sequences
# 2020.06.12
# Shuyu Zheng: shuyu.zheng@helsinki.fi
############################################################

import numpy as np
import pandas as pd
from umi_tools import UMIClusterer
import re
import os
import sys

# check the usage

if len(sys.argv) != 5:
    print("usage: python umitools_cross_sample.py [input_dir] [output_dir] " 
          + "[log_file] [umi_tools_method]")
    sys.exit(1)

## Unpack arguments
input_dir, output_dir, log_file, method = sys.argv[1:]

## Check arguments
if os.path.isdir(input_dir) == False:
    print("Didn't find input directory. \n usage: python umitools_cross_sample.py "
          + "[input_dir] [output_dir] [log_file] [umi_tools_method]")
    sys.exit(2)

methods = ["directional", "adjacency", "percentile", "unique", "cluster"] 

if method not in methods:
    print('The [umi_tools_method] argument is wrong. Available values are: '
          + '"directional", "adjacency", "percentile", "unique", "cluster".'
          + "\n usage: python umitools_cross_sample.py "
          + "[input_dir] [output_dir] [log_file] [umi_tools_method]")
    sys.exit(3)

# Set output, log file path
if not output_dir.endswith("/"):
    output_dir = output_dir + "/"

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
    
# if not os.path.exists(log_file):
#    os.mkdir(log_file)
f = open(log_file + ".txt", 'w+')

def extract_barcode(sequence):
    """
    Extract lineage barcode sequence from reads. The lineage barcode design:
    5’-CTGGGGCACAAGCTTAATTAAGAATTCANNNNTGNNNNACNNNNGANNNNGTNNNNCTAGGGCCTAGAGGGC
    CCGTTTAAAC-3’

    :param sequence: (string) the full sequence for a read.

    :return: (string or None) the subsequence which follow the lineage barcode
            pattern. 
            If there are more than one subsequences are found, all of them are 
            pasted into one string with "," as separator and returned.
            If there is no subsequence detected, this function returns None.
    """
    sequence_search = re.findall('CA.{4}TG.{4}AC.{4}GA.{4}GT.{4}CT', sequence)
    if len(sequence_search)==1:
        sequence_search = sequence_search[0]
    elif len(sequence_search)>=1:
        sequence_search = ','.join(sequence_search)
    else:
        sequence_search = None
    return sequence_search

# Read files from input directory
files = [f for f in os.listdir(input_dir) 
         if os.path.isfile(os.path.join(input_dir, f))]

# Go through the input files and find lineage barcodes from each reads in them.

raw_count = {}
print("Generating naive lineage label count table.")
for fi in files:
    if not os.path.exists(output_dir + fi[:-4]):
        os.mkdir(output_dir + fi[:-4])
    # 1. Read lineage tracing label sequence and counts
    data_path = input_dir + fi
    print(data_path)
    data = pd.read_csv(data_path, sep= ',', header= None)
    data.columns = ['lineage_sequence', 'count','cell_barcode', 'umi']
    f.write('{} has {} reads before removing duplicates, '.format(fi, data.shape))
    data.drop_duplicates(keep = 'first', inplace = True)
    # 2. Extract lineage barcode sequences from raw read sequence
    data['simple_lineage_sequence'] = data['lineage_sequence'].apply(extract_barcode)
    data = data.sort_values(by = ["simple_lineage_sequence", "count"], ascending = False)
    data.drop_duplicates(subset = ['cell_barcode', 'umi'], 
                         keep = 'first', inplace = True)
    raw_count.update({fi: data})
    # data = data.groupby(['cell_barcode', 'lineage_sequence'])['count'].sum().reset_index()
    f.write("{} reads after removing duplicates.\n".format(data.shape))
    data = data.groupby(['cell_barcode', 'simple_lineage_sequence'])['count'].sum().reset_index()
    # data_matrix = data.pivot(index='cell_barcode', columns='simple_lineage_sequence', values='count')
    data.to_csv(output_dir + fi[:-4] + '/naive_sequence_count_df.csv', index = False)
    # data_matrix.to_csv(output_dir + fi[:-4] + '/naive_sequence_count_matrix.csv')

print("Saved naive lineage label count table.")

# Detect UMIs and correct sequence errors (API provided by UMI-tools is used)
# https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/

clusterer = UMIClusterer(cluster_method = method) 
barcode_count_table = pd.concat(list(raw_count.values()))
barcode_count_table = barcode_count_table.groupby(['simple_lineage_sequence']).count().reset_index()
barcode = [el.encode('UTF-8') for el in barcode_count_table['simple_lineage_sequence']]
barcode_table = dict(zip(barcode, barcode_count_table['count']))
clustered_barcode_table = clusterer(barcode_table, threshold = 1)
f.write("Number of lineage barcode clusters: {}\n".format(len(clustered_barcode_table)))

# Generate sequence reference table

tmp = []
for i in range(len(clustered_barcode_table)):
    if len(clustered_barcode_table[i]) == 1:
        represent_sequence = clustered_barcode_table[i][0].decode('UTF-8')
        alternated_sequence_list = []
        alternated_sequence = ""
        total_count = barcode_count_table[barcode_count_table.simple_lineage_sequence == represent_sequence]['count']
        alt_count = ""
        rep_count = total_count
    else:
        alternated_sequence = [x.decode('UTF-8') 
                                  for x in clustered_barcode_table[i]]
        sep_seq_count = barcode_count_table[barcode_count_table\
                                                .simple_lineage_sequence\
                                                .isin(alternated_sequence)]
        sep_count = list(sep_seq_count['count'])
        sep_seq = list(sep_seq_count['simple_lineage_sequence'])
        total_count = sum(sep_count)
        rep_count = max(sep_count)
        max_ind = sep_count.index(rep_count)
        represent_sequence = sep_seq[max_ind]
        alternated_sequence_list = [x for i,x in enumerate(sep_seq) if i != max_ind]
        alternated_sequence = ", ".join(alternated_sequence_list)
        alt_count = [x for i,x in enumerate(sep_count) if i != max_ind]
        alt_count = ", ".join([str(x) for x in alt_count])
    
    tmp.append([represent_sequence, alternated_sequence, rep_count, alt_count, total_count]) 

df = pd.DataFrame(tmp, columns = ['represent_sequence', 'alternated_sequence',
                                  'rep_count', 'alt_count', 'total_count'])

print("Generated sequence reference table.")

df.to_csv(output_dir + "sequence_adjust_reference_table.csv", index = False)

print("Saved sequence reference table.")

# Merge sequences in data and then output count matrix for lineage labels

rep_dic = {i[0].decode('UTF-8') : [x.decode('UTF-8') for x in i]
           for i in clustered_barcode_table}

print("Adjusting sequences")
for fi in files:
    data = raw_count[fi]
    data['adjusted_sequence'] = None
    for i in data.index:
        test_seq = data.simple_lineage_sequence[i]
       # for j in range(len(rep_dic)):
       #     if test_seq == rep_dic[j]:
       #         data.at[i, 'adjusted_sequence'] = test_seq 
       #         break
       #     elif data.simple_lineage_sequence[i] in alt_dic[j]:
       #         data.at[i, 'adjusted_sequence'] = rep_dic[j]
       #         break
        tmp = [k for k, v in rep_dic.items() if test_seq in v]
        if len(tmp) != 0:
            data.at[i, 'adjusted_sequence'] = tmp[0]

    data.to_csv(output_dir + fi[:-4] + '/summary_table.csv', index = False)
    data_adjusted = data.groupby(['adjusted_sequence', 'cell_barcode'])\
                        .count()\
                        .reset_index()
    data_adjusted.drop_duplicates(subset = ['cell_barcode',
                                            'adjusted_sequence'],
                                  keep='first', inplace=True)
    data_adjusted = data_adjusted[['cell_barcode', 'adjusted_sequence', 'count']]
    f.write("Shape of {} UMI count table after merging UMIs from same cluster is {}.\n".format(fi, data_adjusted.shape))
    # data_adjusted_matrix = data_adjusted.pivot(index='cell_barcode', columns='adjusted_sequence', values='count')
    data_adjusted.to_csv(output_dir + fi[:-4] + '/adjusted_sequence_count_df.csv', index = False)
    # data_adjusted_matrix.to_csv(output_dir + fi[:-4] + '/adjusted_sequence_count_matrix.csv')
    print("Saved adjusted count table for file {}.".format(fi))

f.close()
print("Process finished!")
                

