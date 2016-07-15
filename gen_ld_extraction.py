#!/usr/bin/env python
# Copyright (c) 2016 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import logging
import sys
import os
import numpy as np
import pandas as pd
import argparse
from gen_ld_extraction.ld_data import LDData
from gen_ld_extraction.region import Region 

def read_bedfile(bed):
    bed_regions=[]
    with open(bed) as bed:
        for line in bed:
            line = line.split();
            chrom = line[0]
            start = line[1]
            end = line[2]
            name = line[3].strip()
            bed_regions.append(Region(chrom,start,end,name))
    return bed_regions


def read_samples(sam):
    samples = []
    with open(sam) as samf:
        for i,line in enumerate(samf):
            if i > 1:
                samples.append(line.split()[1])
    return samples


def read_keep(keep_f, samples):
    keep=set()
    with open(keep_f,'r') as keepf:
        for line in keepf.readlines():
            #indiv = '_'.join(line.strip().split()[i] for i in (0,1))
            indiv = line.strip().split()[0]
            keep.add(indiv)
            #print keep
    keep = [(i) for i, e in enumerate(samples) if e in keep]
    return keep 


#def read_data():
#    snpids = []
#    A1s = []
#    A2s = []
#    data = pd.read_csv(sys.stdin, sep=" ", header=None)
#
#    snp=0
#    snpids = []
#    A1s = []
#    A2s = []
#    for line in sys.stdin:
#        snpids.append(raw[0])
#        A1s.append(raw[1])
#        A2s.append(raw[2])
#        for i in xrange(nsamples):
#            ii = 2*i+3
#            data[snp,i] = float(raw[ii])*2.0+float(raw[ii+1])
#        snp+=1
#
#
#def sample_subset_data():
#    keepis = [i for i, e in enumerate(samples) if e in keep or i in range(6)]
#    nkeepsamples = len(keepis)
#    return data[:, keepis]
#
#
#def recalc_mafs(data):
#    return None
#
#
#def snps_subset_data(data, maf=0):
#    data.ix[:,5] = data.ix[:,6:].apply(np.mean, axis=1)/2
#    data = data.ix[(data.ix[:,5]>0.01) & (data.ix[:,5]<0.99)]
#
#
#
#def calc_ldmatrix():
#    x = np.empty([data.shape[0],data.shape[0]])
#    x = np.corrcoef(data.ix[:,6:])
#
#
#def write_ldmatrix(OUT, matrix):
#    with open(OUT+".snps", 'w') as snpsout:
#
#    #with open(OUT+".ld", 'w') as ldout
#        np.savetxt(fout, x)



def main():
    parser = argparse.ArgumentParser(description="CMC LD matrix extraction")
    parser.add_argument("-g", "--gen_fileroot", dest="gen_fileroot", required=True)
    parser.add_argument("-s", "--samples", dest="samplefile_for_gen_files", required=True)
    parser.add_argument("-k", "--keep", dest="keep_file", required=True)
    parser.add_argument("-o", "--out", dest="output_fileroot", required=True)
    parser.add_argument("-b", "--bed_intervals", dest="bed_file", required=True)
    parser.add_argument("-m", "--maf", dest="maf_filter", default=0.01)
    args = parser.parse_args()
    args.maf_filter  = float(args.maf_filter)

    try:
        os.mkdir(args.output_fileroot)
    except OSError:
        pass
    logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    
    fileHandler = logging.FileHandler(os.path.join(args.output_fileroot,args.bed_file+".log"),mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    logging.info("Welcome to CMC LD matrix extraction")

    logging.info("Starting to load Bed regions")
    bed_regions = read_bedfile(args.bed_file)
    logging.info("Loaded {0} bed regions".format(len(bed_regions)))

    logging.info("Starting to load sample file") 
    samples = read_samples(args.samplefile_for_gen_files)
    logging.info("Loaded {0} samples".format(len(samples)))
    
    logging.info("Starting to load keep file") 
    keep_samples = read_keep(args.keep_file, samples)
    logging.info("Loaded {0} samples to keep".format(len(keep_samples)))
    
    for region in bed_regions:
        ld_data = LDData(region, args.gen_fileroot,  keep_samples)
        logging.info("Loading gene files(s) and converting to dosages")
        try:
            ld_data.load_gen_and_generate_dosages()
            logging.info("Loaded gen file(s) and converted to dosages format")
            ld_data.filter_maf(args.maf_filter)
            ld_data.calc_ld()
            ld_data.write_outputs(args.output_fileroot)
            logging.info("Completed processing {0}".format(str(region)))
        except Warning as e:
            logging.error(e)
            pass

if __name__=="__main__":
    main()
