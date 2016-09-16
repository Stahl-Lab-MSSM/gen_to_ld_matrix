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


import glob
import logging
import os
import pandas as pd
import numpy as np
import sys

class GenData(object):
    """
        Dosage data object.
    """
    def __init__(self,reg,genage_file_root, keep):
        self._start = reg.start
        self._end = reg.end
        self._chrom = reg.chrom
        self._keep_list = keep
        self._file_root = genage_file_root
        self._gen = None
    
    @property
    def start(self):
        return self._start
    @property 
    def end(self):
        return self._end
    @property 
    def chrom(self):
        return self._chrom
    @property
    def keep_list(self):
        return self._keep
    @property
    def gen(self):
        return self._gen
    @property
    def file_root(self):
        return self._file_root
    @property
    def keep_list(self):
        return self._keep_list
    def update_gen(self, gen):
        self._gen = gen 

    def get_dosage_frame(self):
        logging.info("Starting to create the dosages frame")
        number_of_samples = (len(self.gen.columns) -5 )/3
        dosage_frame = pd.DataFrame(index=self.gen.index, columns=range(number_of_samples +7))
        for i in range(5):
            dosage_frame.iloc[:,i]=self.gen.iloc[:,i]
        af = 0.0 
        for j in range(number_of_samples):
            if(j % 10 == 0):
                logging.info("Processed {0} samples".format(j))
            jj = 3* j + 5
            dosage_frame.iloc[:,j+7] = self.gen.iloc[:,jj] * 2.0 + self.gen.iloc[:,jj+1]
        dosage_frame.iloc[:,5] = dosage_frame.iloc[:,7:].apply(np.mean,axis=1) /2
        dosage_frame.iloc[:,6] = dosage_frame.iloc[:,7:].apply(np.var,axis=1) 
        logging.info("Calculated the variance")
        logging.info("Created the dosages frame")
        return dosage_frame

    def sample_filter(self):
        logging.info("Removing samples not in the keep file")
        #sample_filter = [False] * self.gen
        keep_columns = range(5)
        # Keep list zero indexed and 3 values per sample so .....
        # i * 3 +5 = first sample genotype
        for i in self.keep_list:
            keep_columns.append(i*3 + 5 )
            keep_columns.append(i*3 + 1 + 5 )
            keep_columns.append(i*3 + 2 + 5 )
        self.update_gen(self.gen.ix[:,keep_columns])
        logging.info("Dosage matrix after removing samples, shape = {0}".format(self.gen.shape))


    def load_gen(self, gen_file=None):
        logging.info("Attempting to load region = chr{0}:{1}-{2}".format(self.chrom, self.start,self.end))
        start_mb = self.start/1e6
        end_mb= self.end/1e6

        if(gen_file is None):
            logging.info("File root {0}".format(self.file_root))
            gen_path = os.path.abspath(os.path.dirname(self.file_root))
            gen_basename = os.path.basename(self.file_root)
            cmc_files = glob.glob(os.path.join(gen_path, gen_basename + "{0}_*.gen".format(self.chrom)))
            try:
                cmc_files[0]
            except KeyError:
                logging.info("Expecting gene files in File root {0} to be gzipped (*.gen.gz)".format(self.file_root))
                cmc_files = glob.glob(os.path.join(gen_path, gen_basename + "{0}_*gen.gz".format(self.chrom)))
                pass
            cmc_files = sorted(cmc_files, key= lambda x: int(x.split('.')[1].split('-')[0]))
            files_to_load = []
        else:
            cmc_files = [gen_file]
        for cmc_file in cmc_files:
            cmc_start = int(cmc_file.split('.')[1].split('-')[0])
            cmc_end = int(cmc_file.split('.')[1].split('-')[1].replace("Mb",""))
            if (cmc_start <= start_mb and cmc_end >= start_mb) or (cmc_start <= end_mb and cmc_end >= start_mb) : 
                files_to_load.append(cmc_file)
            if cmc_end >= end_mb:
                break
        logging.info("Loading gen_files {0}".format(" ".join(files_to_load)))
        if len(files_to_load) < 1:
            #logging.error("Could not find any gen files for chr{0}:{1}-{2}".format(self.chrom, self.start, self.end))
            #sys.exit(1)
            raise Warning("Could not find any gen files for chr{0}:{1}-{2}".format(self.chrom, self.start, self.end))

        for i, gen_file in enumerate(files_to_load):
            logging.info("Loading {0}".format(gen_file))
            if i < 1: 
                data = pd.read_csv(gen_file, sep=" ", header=None)
                logging.info("Gen matrix before filtering SNPs: shape = {0}".format(str(data.shape)))
                data= data[(data.ix[:,2] >= self.start) & (data.ix[:,2] <= self.end)]
                logging.info("Gen matrix after filtering SNPs: shape = {0}".format(str(data.shape)))
            else:
                data_tmp  = pd.read_csv(gen_file, sep=" ", header=None)
                data_tmp = data_tmp[(data_tmp.ix[:,2] >= self.start) & (data_tmp.ix[:,2] <= self.end)]
                data = pd.concat([data,data_tmp])
        # convert gene_data to dos_data
        self.update_gen(data)
        #data.to_csv(os.path.join(gen_file), sep=" ", header=False, index=False)
