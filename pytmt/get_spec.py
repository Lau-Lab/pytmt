# -*- coding: utf-8 -*-

""" Reads in mzml file using pymzml and get list of ms2 scans """

import logging
import pymzml as mz


class Mzml(object):
    """ Mzml class. """

    def __init__(
            self,
            path: str,
            precision: int = 20,
            logger: logging.Logger = None,
    ) -> None:
        """
        This class reads mzml files using pymzml

        :param path: path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        :param precision: Int determines precision of reading as well as mass tolerance of peak integration (ppm)
        """

        self.path = path    # path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        self.ms2data = {}    # dictionary of ms2 spectra
        self.ms3data = {}    # dictionary of ms3 spectra
        self.prec_idx = {}  # dictionary of precursor scan id
        self.rt_idx = {}    # dictionary of retention times
        self.mslvl_idx = {} # dictionary of ms levels
        self.precision = precision  # integer determines precision of reading as well as mass tolerance of peak integration (ppm)
        self.logger = logger if logger else logging.getLogger(__name__)     # logger

    def parse_mzml_ms2(self) -> None:
        """
        Read the mzml file and create data dictionary for all ms2 peaks
        :return:
        """

        run = mz.run.Reader(self.path)

        # 2023-10-04: added ms3 support
        run.ms_precisions = {1: self.precision * 1e-6,
                             2: self.precision * 1e-6,
                             3: self.precision * 1e-6}


        for n, spec in enumerate(run):

            self.mslvl_idx[n + 1] = spec.ms_level
            self.rt_idx[n + 1] = spec.scan_time

            if spec.ms_level == 2:
                self.ms2data[n + 1] = spec.peaks("centroided")
                self.prec_idx[n + 1] = spec.selected_precursors[0].get('precursor id')

            elif spec.ms_level == 3:
                self.ms3data[n + 1] = spec.peaks("centroided")
                self.prec_idx[n + 1] = spec.selected_precursors[0].get('precursor id')

            else:
                continue

        self.logger.info(f'Parsed {n + 1} spectra from file {self.path}')

        return None
