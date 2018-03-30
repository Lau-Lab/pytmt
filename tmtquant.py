import pymzml as mz
import os.path

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


