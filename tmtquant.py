import pymzml as mz
import os.path
import re

#import os.path
#import pandas as pd
#import scipy.integrate
from time import time


class Mzml(object):
    """
    This class handles the Mzml files to be read
    """

    def __init__(self, path):
        """

        :param path:
        """

        self.path = os.path.join(path)
        self.msdata = mz.run.Reader(self.path, MS1_Precision=20e-6, MSn_Precision=20e-6)

    def get_reporter_from_scan_id(self, reporter_list, spectrum_id):
        """

        :param reporter_list:
        :param spectrum_id:
        :return:
        """
        pass




mzid_loc = "./test_mzml_1"
perc_loc = "./test_perc_1"

assert os.path.isdir(perc_loc)

# List all files in the percolator directory ending with target.psms.txt
perc_psms_file = [f for f in os.listdir(perc_loc) if f.endswith('target.psms.txt')]

assert len(perc_psms_file) == 1, 'Check percolator output directory has 1 target.psms.txt'

