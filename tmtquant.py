import pymzml as mz
import os.path

import re

#import os.path
import pandas as pd
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




mzml_loc = "./test_mzml_1"
perc_loc = "./test_perc_1"

assert os.path.isdir(perc_loc)

# List all files in the percolator directory ending with target.psms.txt
perc_psms_files = [f for f in os.listdir(perc_loc) if f.endswith('target.psms.txt')]

assert len(perc_psms_files) == 1, 'Check percolator output directory has 1 target.psms.txt'

# Read the psms file.
psms_df = pd.read_csv(filepath_or_buffer=os.path.join(perc_loc, perc_psms_files[0]),
                      sep="\t")

# Get all the file indices in the run
file_indices = list(set(psms_df['file_idx']))

# Check that the number of mzMLs in the mzml folder is the same as the maximum of the ID's file_idx column
# Note this will throw an error if not every fraction results in at least some ID, but we are going to
# ignore this for now

mzml_files = [f for f in os.listdir(mzml_loc) if re.match('^.*.mzML', f)]
# Sort by names
mzml_files.sort()

# Throw an error if there is no mzML file in the mzml directory
if len(mzml_files) == 0:
    try:
        raise Exception('No mzML files in the specified directory')
    except ValueError as e:
        print('Unknown error')

assert len(mzml_files) == max(file_indices) + 1, 'Number of files not matching'

# Loop through each file index, then

for idx in file_indices:

    # Make a sub_df with current index being considered
    sub_df = psms_df[psms_df['file_idx'] == idx]

    # Arrange by scan number
    sub_df = sub_df.sort_values(by='scan')

    # Filter by q value
    sub_df_filter = sub_df[sub_df['percolator q-value'] <= 0.05]
    sub_df_filter = sub_df_filter.reset_index()

    # Open the mzML file
    msdata = mz.run.Reader(path=os.path.join(mzml_loc, mzml_files[idx]),
                  MS1_Precision=20e-6, MSn_Precision=20e-6)


    # Loop through each qualifying row in sub_df_filtered

    for i in range(len(sub_df_filter)):

        # scan number
        scan_num = sub_df_filter.ix[i, 'scan']

        # Get the spectrum based on the spectrum number
        try:
            spectrum = msdata[scan_num]

        except KeyError:
            print("Spectrum not found")

        tmt_list = [126.127726,
                    127.124761,
                    127.131081,
                    128.128116,
                    128.134436,
                    129.131471,
                    129.137790,
                    130.134825,
                    130.141145,
                    131.138180
                    ]

        tmt_intensities = []

        for reporter in tmt_list:

            matchList = spectrum.has_peak(reporter)

            if matchList:
                for mz, I in matchList:
                    tmt_intensities.append([scan_num, i, I, mz])
            else:
                tmt_intensities.append([scan_num, i, 0, 0])


