# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 15:08:41 2020

@author: cyclopenta
"""


import pandas as pd
import os
import re
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description='整合各样本的STAR比对信息（默认为./log_final_out目录下所有*Log.final.out文件）: \
                                              python inf_star_map_ver0.1.py \
                                              -i input dir \
                                              -o output')

parser.add_argument('-i', '--input_dir', type=str, default='./log_final_out', help='url to /log_final_out/')
parser.add_argument('-o', '--output_prefix', type=str, required=True, help='output_file_prefix')

args = parser.parse_args()

#all the items listed here
'''
['Started job on ', '\tAug 07 13:59:41\n']
['Started mapping on ', '\tAug 07 14:07:16\n']
['Finished on ', '\tAug 07 14:25:01\n']
['Mapping speed, Million of reads per hour ', '\t169.78\n']
['']
['Number of input reads ', '\t50227692\n']
['Average input read length ', '\t150\n']
['UNIQUE READS:\n']
['Uniquely mapped reads number ', '\t1545903\n']
['Uniquely mapped reads % ', '\t3.08%\n']
['Average mapped length ', '\t140.03\n']
['Number of splices: Total ', '\t157355\n']
['Number of splices: Annotated (sjdb) ', '\t107944\n']
['Number of splices: GT/AG ', '\t116913\n']
['Number of splices: GC/AG ', '\t2546\n']
['Number of splices: AT/AC ', '\t240\n']
['Number of splices: Non-canonical ', '\t37656\n']
['Mismatch rate per base, % ', '\t1.39%\n']
['Deletion rate per base ', '\t0.13%\n']
['Deletion average length ', '\t3.11\n']
['Insertion rate per base ', '\t0.07%\n']
['Insertion average length ', '\t1.17\n']
['MULTI-MAPPING READS:\n']
['Number of reads mapped to multiple loci ', '\t3954028\n']
['% of reads mapped to multiple loci ', '\t7.87%\n']
['Number of reads mapped to too many loci ', '\t3810\n']
['% of reads mapped to too many loci ', '\t0.01%\n']
['UNMAPPED READS:\n']
['Number of reads unmapped: too many mismatches ', '\t0\n']
['% of reads unmapped: too many mismatches ', '\t0.00%\n']
['Number of reads unmapped: too short ', '\t44719792\n']
['% of reads unmapped: too short ', '\t89.03%\n']
['Number of reads unmapped: other ', '\t4159\n']
['% of reads unmapped: other ', '\t0.01%\n']
['CHIMERIC READS:\n']
['Number of chimeric reads ', '\t0\n']
['% of chimeric reads ', '\t0.00%\n']
'''

# initialization
total_info = Counter()
interest_info = Counter()
temp_dic = Counter()
#define what you want
interest_ini = ['Filename', 'Sample', 'Uniquely mapped reads %', '% of reads mapped to multiple loci', '% of reads unmapped: too short', 'Number of input reads']
for ini in interest_ini:
    interest_info[ini] ={}
    
#open dir and read files    
file_path = args.input_dir
file_total = os.listdir(file_path)
file_total = sorted(file_total)

for f in file_total:
    #print('f',f)
    if not os.path.isdir(f):
        file = open(file_path+'/'+f)
        onefile_info = Counter()
        onefile_info['Filename'] = re.sub('Log.final.out', '', f)
        onefile_info['Sample'] = re.sub('Log.final.out', '', f).split('_')[0]
        for line in file:
            line = line.lstrip()
            i = line.split('|')
            # prevent list index error
            if len(i) > 1:
                onefile_info[i[0].rstrip()] = i[1].lstrip('\t').rstrip()
        #print(onefile_info)
        total_info[f] = onefile_info
#print(total_info)
    
for k,v in total_info.items():
    for m,n in interest_info.items():
        interest_info[m][k]=v[m]
# print(interest_info)

final_df = pd.DataFrame(interest_info)
final_df.columns = ['Filename','Sample','Uniquely mapped reads %','Reads mapped to multiple loci %','% of reads unmapped: too short','Number of input reads']
final_df['Mapped reads %'] = final_df['Uniquely mapped reads %'].str.strip('%').astype(float)+final_df['Reads mapped to multiple loci %'].str.strip('%').astype(float)
final_df['Mapped reads %'] = final_df['Mapped reads %'].apply(lambda x: format(x, '.2f')).astype(str)+'%'
final_df = final_df[['Filename','Sample','Uniquely mapped reads %','Reads mapped to multiple loci %','Mapped reads %','% of reads unmapped: too short','Number of input reads']]
# output_name = input('Enter the output name (date) :\n')
final_df.to_csv(args.output_prefix+'_RNA_MAP.csv', index=False, sep=',')

