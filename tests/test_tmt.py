# -*- coding: utf-8 -*-

""" Tests """

import unittest
import os
import ftplib
import pandas as pd
from tqdm import tqdm


class MzmlTest(unittest.TestCase):
    """
    Test cases involving reading Mzml files
    """

    def setUp(self):
        """
        Check that files exist, otherwise download from ProteomeXchange

        :return:
        """

        mzml_files = ['20180416_StemCell_TMT_Block6_F7.mzML', '20180416_StemCell_TMT_Block6_F8.mzML']

        test_path = os.path.join('tests', 'data', 'mzml')

        file1_exists = os.path.isfile(os.path.join(test_path, mzml_files[0]))
        file2_exists = os.path.isfile(os.path.join(test_path, mzml_files[1]))

        if file1_exists and file2_exists:
            return True

        if not os.path.exists(test_path):
            os.makedirs(test_path)

        ftp = ftplib.FTP("ftp.pride.ebi.ac.uk")
        ftp.login(user='', passwd='')
        ftp.cwd('/pride/data/archive/2019/12/PXD013426/')

        for px_file in mzml_files:
            size = ftp.size(px_file)
            host_file = os.path.join(test_path, px_file)

            try:
                with open(host_file, 'wb') as local_file:

                    with tqdm(total=size,
                              unit='B', unit_scale=True, unit_divisor=1024,
                              ) as pbar:
                        pbar.set_description("Downloading test data from ProteomeXchange")

                        def callback(data):
                            pbar.update(len(data))
                            local_file.write(data)

                        ftp.retrbinary('RETR {}'.format(px_file), callback)

            except ftplib.error_perm:
                print('error FTP')
                pass

        ftp.quit()

    def tearDown(self):
        pass

    def test_that_percolator_file_exists(self):
        """
        Check that there is one percolator file
        """

        percolator_path = os.path.join('tests', 'data', 'percolator')
        id_files = [f for f in os.listdir(percolator_path) if f.endswith('target.psms.txt')]

        self.assertEqual(len(id_files), 1)

    def test_that_percolator_file_opens(self):
        """
        Check that percolator file can be opened by pandas and has two Mzml indices
        """

        percolator_path = os.path.join('tests', 'data', 'percolator')
        id_files = [f for f in os.listdir(percolator_path) if f.endswith('target.psms.txt')]

        id_df = pd.read_csv(filepath_or_buffer=os.path.join(percolator_path, id_files[0]),
                            sep='\t')

        # Get all the file indices in the Percolator results file.
        file_indices = list(set(id_df['file_idx']))

        self.assertEqual(len(file_indices), 2)


