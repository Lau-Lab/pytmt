# -*- coding: utf-8 -*-


""" Reads in mzml file using pymzml and get list of ms2 scans """

import logging
import os
import pymzml as mz


class Mzml(object):

    def __init__(
            self,
            path,
            precision
    ):
        """
        This class reads mzml files using pymzml

        :param path: path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        :param precision: integer determines precision of reading as well as mass tolerance of peak integration (ppm)
        """

        self.path = path
        self.msdata = {}
        self.rt_idx = {}
        self.mslvl_idx = {}
        self.precision = precision

        self.logger = logging.getLogger('pytmt.mzml')

    def parse_mzml_ms2(self):
        """
        Read the mzml file and create data dictionary for all ms2 peaks
        :return:
        """

        # 2022-03-28 try to open either mzML or mzML.gz
        if os.path.exists(self.path + '.mzML'):
            open_path = self.path + '.mzML'
        elif os.path.exists(self.path + '.mzML.gz'):
            open_path = self.path + '.mzML.gz'
        else:
            raise FileNotFoundError

        run = mz.run.Reader(open_path,
                            MS_precision={
                                1: self.precision*1e-6,
                                2: self.precision*1e-6
                            })



        for n, spec in enumerate(run):

            #if n % 1000 == 0:
            #    print(
            #        'Loading spectrum {0} at retention time {scan_time:1.2f}'.format(
            #           spec.ID,
            #            scan_time=spec.scan_time
            #        )
            #   )

            self.mslvl_idx[n + 1] = spec.ms_level
            self.rt_idx[n + 1] = spec.scan_time

            if spec.ms_level == 2:
                self.msdata[n + 1] = spec.peaks("centroided")

        self.logger.info(
                'Parsed {0} spectra from file {1}'.format(
                    n + 1,
                    self.path)
                )

        return True

