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
        self.msdata = {}    # dictionary of ms2 spectra
        self.rt_idx = {}    # dictionary of retention times
        self.mslvl_idx = {} # dictionary of ms levels
        self.precision = precision  # integer determines precision of reading as well as mass tolerance of peak integration (ppm)
        self.logger = logger if logger else logging.getLogger(__name__)     # logger

    def parse_mzml_ms2(self) -> None:
        """
        Read the mzml file and create data dictionary for all ms2 peaks
        :return:
        """

        run = mz.run.Reader(self.path,
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

        self.logger.info(f'Parsed {n + 1} spectra from file {self.path}')

        return None

