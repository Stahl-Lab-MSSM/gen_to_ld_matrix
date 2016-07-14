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
        self._gen_data = None
    
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
    def dosages(self):
        return self._gen_data
    @property
    def file_root(self):
        return self._file_root
    @property
    def keep_list(self):
        return self._keep_list
    def update_gen(self, dosages):
        self._gen_data = dosages 

    def get_dosage_frame(self):
        logging.info("Starting to create the dosages frame")
        number_of_samples = (len(self.dosages.columns) -5 )/3
        dosage_frame = pd.DataFrame(index=self.dosages.index, columns=range(number_of_samples +6))
        for i in range(len(self.dosages.index)):
            if( i % 1000 == 0):
                logging.info("Processed {0} rows".format(i))
            meta_info = []
            for j in range(5):
                meta_info.append(self.dosages.iat[i,j]) 
            af = 0.0 
            offset = 5
            dosages = []
            for j in range(number_of_samples):
                jj = 3* j + 5  
                tmp_dos = self.dosages.iat[i,jj] * 2.0 + self.dosages.iat[i,jj+1]
                dosages.append(tmp_dos)
                af += tmp_dos
            af = af/(number_of_samples * 2)
            if (af > 0.5):
                af = 1 - af
            meta_info.append(af)
            meta_info.extend(dosages)
            dosage_frame.iloc[i] = meta_info
        logging.info("Created the dosages frame")
        return dosage_frame

    def sample_filter(self):
        logging.info("Removing samples not in the keep file")
        sample_filter = [False] * self.dosages
        keep_columns = range(5)
        # Keep list zero indexed and 3 values per sample so .....
        # i * 3 +5 = first sample genotype
        for i in self.keep_list:
            keep_columns.append(i*3 + 5)
            keep_columns.append(i*3 + 1 + 5)
            keep_columns.append(i*3 + 2+ 5 )
        self.update_gen(self.dosages.ix[:,keep_columns])
        logging.info("Dosage matrix after removing samples, shape = {0}".format(self.dosages.shape))


    def load_gen(self):
        logging.info("Attempting to load region = chr{0}:{1}-{2}".format(self.chrom, self.start,self.end))
        start_mb = (self.start/1e6)//5* 5
        end_mb= (self.end/1e6+5)//5 * 5
        logging.info("File root {0}".format(self.file_root))
        cmc_files = glob.glob(os.path.join(os.path.abspath(self.file_root), "CM6-pos-chr{0}*gen".format(self.chrom))) 
        cmc_files = sorted(cmc_files, key= lambda x: int(x.split('.')[1].split('-')[0])) 
        files_to_load = []
        for cmc_file in cmc_files:
            cmc_start = int(cmc_file.split('.')[1].split('-')[0])
            cmc_end = int(cmc_file.split('.')[1].split('-')[1].replace("Mb",""))
            if cmc_end > end_mb:
                break
            if cmc_start >= start_mb: 
                files_to_load.append(cmc_file)
        logging.info("Loading gen_files {0}".format(" ".join(files_to_load)))
        if len(files_to_load) < 1:
            logging.error("Could not find any gen files for chr{0}:{1}-{2}".format(self.chrom, self.start, self.end))
            sys.exit(1)

        for i, gen_file in enumerate(files_to_load):
            logging.info("Loading {0}".format(gen_file))
            if i < 1: 
                data = pd.read_csv(gen_file, sep=" ", header=None)
                logging.info("Dosage matrix before filtering SNPs: shape = {0}".format(str(data.shape)))
                data= data[(data.ix[:,2] >= self.start) & (data.ix[:,2] <= self.end)]
                logging.info("Dosage matrix after filtering SNPs: shape = {0}".format(str(data.shape)))
            else:
                data_tmp  = pd.read_csv(gen_file, sep=" ", header=None)
                data_tmp = data[(data.ix[:,2] >= self.start) & (data.ix[:,2] <= self.end)]
                data = pd.concat(data_tmp)
        self.update_gen(data)
        #data.to_csv(os.path.join(gen_file), sep=" ", header=False, index=False)
