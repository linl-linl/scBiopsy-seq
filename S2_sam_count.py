# ！usr/bin/env python
# -*- coding:utf-8 -*-

# 用htseq-count，由Aligned.out.sam计算read counts
# 调整read counts文件为count表达

import os
import re
import argparse

parser = argparse.ArgumentParser(description='由目录下的(指定的)Aligned.out.sam提取read counts，生成表达矩阵: \
                                              python cluster_sam_count.py \
                                              -s list_sam \
                                              -g url_gtf')

parser.add_argument('-s', '--sam', type=str, nargs='+', help='a sam file(s) with suffix".sam" or a list of sam')
parser.add_argument('-g', '--gtf', type=str, default='/home/disk/linl/reference/GRCh38/Homo_sapiens.GRCh38.105.gtf', help='gtf file for htseq-count')

args = parser.parse_args()


# save sam_url to list_sam
if args.sam:
    list_sam = []
    for ss in args.sam:
        if ss.endswith('.sam'):
            list_sam.append(ss)
        else:
            file = open(ss, 'r')
            for i in file:
                list_sam.append(i.rstrip())
            file.close()
    list_sam = sorted(list(set(list_sam)))
else:
    os.system("ls *Aligned.out.sam > list_sam_cnt.tmp.txt")
    #os.system("ls *depth.sam > list_sam_htseq.tmp.txt")
    #os.system("ls *.sam > list_sam_htseq.tmp.txt")
    list_sam = []
    file_sam = open('list_sam_cnt.tmp.txt', 'r')
    for i in file_sam:
        list_sam.append(i.rstrip())
    os.system("rm list_sam_cnt.tmp.txt")

# htseq-count: output "gene\tcount"
gtf = args.gtf
list_counts = []
for sam in list_sam:
    prefix = re.sub('Aligned.out.sam', '', sam)
    prefix = re.sub('depth.sam', '', prefix)
    prefix = re.sub('.sam', '', prefix)
    htseq_out1 = prefix+'_union_counts.txt'
    htseq_out2 = prefix+'_nonempty_counts.txt'
    list_counts.append(htseq_out1)
    list_counts.append(htseq_out2)
    # htseq-count
    os.system("htseq-count -f sam -r name -s no -a 10 -t exon -i gene_id -m union %s %s > %s" % (sam, gtf, htseq_out1))
    os.system("htseq-count -f sam -r name -s no -a 10 -t gene -i gene_id -m intersection-nonempty %s %s > %s" % (sam, gtf, htseq_out2))

# rm statistical information from htseq_out
# add colname to htseq_out
# output single cell dge
## CountsToCount
os.system("ls *_counts.txt > list_cnts_cnt.tmp.txt")
list_cnts = []
file_cnts = open('list_cnts_cnt.tmp.txt', 'r')
for i in file_cnts:
    list_cnts.append(i.rstrip())
os.system("rm list_cnts_cnt.tmp.txt")
## list_counts = list_cnts
for cs in list_counts:
    prefix = re.sub('_counts.txt', '', cs)
    count_filename = prefix+'_count.txt'

    file_counts = open(cs, 'r')
    file_count = open(count_filename, 'w')
    counts_cont = file_counts.readlines()
    count_cont = counts_cont[:-5]
    file_count.write('GENE\t'+prefix.split('_')[0]+'\n')
    file_count.write(''.join(count_cont))
    file_counts.close()
    file_count.close()
    print(count_filename+" finished!!!")

print("All Finished!!!")

