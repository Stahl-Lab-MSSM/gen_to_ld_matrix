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

import os
from gen_ld_extraction.gen_data import GenData 
import numpy as np
import logging

class LDData(object):
    """
       Highest level LDData object. 
    """

    def __init__(self, reg, gen_file_root, keep):
        self._reg = reg
        self._gen_file_root = gen_file_root
        self._keep = keep

    @property
    def reg(self):
        return self._reg
    @property
    def gen_file_root(self):
        return self._gen_file_root
    @property
    def samples(self):
        return self._samples
    @property
    def keep(self):
        return self._keep
#    @property
#    def gen_data(self):
#        return self._gen_data
    @property
    def dosages(self):
        return self._dosages
    @property
    def ld_matrix(self):
        return self._ld_matrix
    def update_ld_matrix(self,ld_matrix):
        self._ld_matrix = ld_matrix

    def update_dosages(self, dosages):
        self._dosages = dosages 

    def filter_maf(self,maf):
        logging.info("Filtering by MAF {0} ".format(maf))
        self.update_dosages(self.dosages.ix[(self.dosages.ix[:,5]>maf) & (self.dosages.ix[:,5]<(1-maf))])
        logging.info("Shape after MAF filtering {0}".format(self.dosages.shape))

    def calc_ld(self):
        logging.info("Calculating LD matrix")
        self.dosages.to_csv("Test2.txt", sep=" ", header=None, index=False)
        ld_matrix=np.corrcoef(self.dosages.ix[:,7:].astype(float).values)
        logging.info("Calculated LD matrix") 
        self.update_ld_matrix(ld_matrix)

    def load_gen_and_generate_dosages(self, gene_file=None):
        # initialise dosages
        gen_data = (GenData(self.reg, self.gen_file_root, self.keep, gen_file=None))
        gen_data.load_gen()
        gen_data.sample_filter()
        self.update_dosages(gen_data.get_dosage_frame())


    def write_outputs(self, output_root):
        meta_frame = self.dosages.ix[:,1:6]
        meta_frame.columns = ["RSID","POS","A1","A2","FREQ1","GVAR"]
        meta_frame.to_csv(os.path.join(output_root, self.reg.name + '.snpdat'),index=False, sep= " ")
        np.savetxt(os.path.join(output_root, self.reg.name+".ld"),self.ld_matrix)
        self.dosages.to_csv(os.path.join(output_root, self.reg.name + ".dos"), index=False, sep=" " )




