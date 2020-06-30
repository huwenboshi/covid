#!/usr/bin/python
# (c) 2018-2023 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys, logging, time
from src.estimation import *

# main function
def main():

    # get command line argument and initialize log
    args = get_command_line()

    # parse command line 
    argmap = check_command_line(args)

    # load frequency file
    freq = load_frqfile(args.frqfile, start_chrom, stop_chrom)

    # load seg annotations
    annot_snps,all_seg_annot,seg_annot_mat = load_annot(argmap['seg-annot'],
        start_chrom, stop_chrom)
    annot_idx = intersect_snp_estimation(annot_snps['SNP'], freq, 0.05)

def load_annot(prefix_list, start_chrom, end_chrom,
    annot_start_idx=4, snp_idx=2):
    """
    Load annotation of all chromosomes
    """

    # load all snps snps
    all_snp = []
    for i in xrange(start_chrom, end_chrom+1):
        prefix = prefix_list[0]
        filename = '{}{}.annot.gz'.format(prefix, i)
        tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
            na_filter=False, memory_map=True, usecols=['SNP'])
        all_snp.append(tbl)
    all_snp = pd.concat(all_snp, axis=0, ignore_index=True)
    tot_nsnp = all_snp.shape[0]

    # get annotations
    all_annot = []
    annot_list = []
    for i in xrange(len(prefix_list)):
        prefix = prefix_list[i]
        filename = '{}{}.annot.gz'.format(prefix, start_chrom)
        with gzip.open(filename) as f:
            line = f.readline().strip()
            tmp = line.split()[annot_start_idx:]
            all_annot += tmp
            annot_list.append(tmp)
    tot_nannot = len(all_annot)
    
    # load the annot into memory
    all_annot_mat = np.zeros((tot_nsnp, tot_nannot), dtype=np.float32)
    idx_c = 0
    for k in xrange(len(prefix_list)):
        prefix = prefix_list[k]
        idx_r = 0
        for i in xrange(start_chrom, end_chrom+1):
            filename = '{}{}.annot.gz'.format(prefix, i)
            tmp = annot_list[k]
            dt_load = dict(zip(tmp, [np.float32]*len(tmp)))
            tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
                na_filter=False, memory_map=True, usecols=tmp, dtype=dt_load)
            nrow, ncol = tbl.shape
            all_annot_mat[idx_r:idx_r+nrow,idx_c:idx_c+ncol] = tbl[tmp].values
            idx_r += nrow
        idx_c += ncol
    
    # convert annot to data frame
    all_annot = pd.DataFrame({'ANNOT': all_annot})

    return all_snp, all_annot, all_annot_mat

# intersect snp
def intersect_snp_estimation(annot_snps, freq, min_maf):
    annot_idx = np.where(freq['MAF']>min_maf)[0]
    return annot_idx

# load minor allele frequency file
def load_frqfile(frqfile_fnm, start_chrom, stop_chrom):
    
    all_frq = []
    for i in xrange(start_chrom, stop_chrom+1):
        all_frq.append(pd.read_table('{}{}.frq'.format(frqfile_fnm,i),
            delim_whitespace=True))
    
    all_frq = pd.concat(all_frq, axis=0, ignore_index=True)

    return all_frq

# create blocks for block jackknife
def create_block(start_idx, stop_idx, nblock):

    block_size = int(np.ceil(float(stop_idx-start_idx+1)/float(nblock)))
    cuts = range(start_idx, stop_idx, block_size)
    blocks = []
    
    for i in range(len(cuts)-1):
        start = cuts[i]
        stop = cuts[i+1]
        blocks.append(np.arange(start, stop))
        
    blocks.append(np.arange(cuts[len(cuts)-1], stop_idx+1))

    return blocks

# check command line
def check_command_line(args):

    # parse command line
    argmap = dict()
    for arg in vars(args):
        argmap[arg] = getattr(args, arg)

    return argmap

# get command line input
def get_command_line():
    
    # create the help document
    parser = argparse.ArgumentParser(description='Get correlation between'
        'chisq statistics and SEG annotation')

    parser.add_argument('--seg-annot', dest='seg-annot', type=str,
        required=False, nargs='+', help='Baseline annotation file')

    parser.add_argument('--use-chrom', dest='use-chrom', type=int,
        required=False, default=[1, 22], nargs=2,
        help='Specific which chromosomes to use')

    parser.add_argument('--frqfile', dest='frqfile', type=str, nargs=2,
        required=False, help='Prefix of the frequency files')

    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

# execute the main function
if(__name__ == '__main__'):
    main()
