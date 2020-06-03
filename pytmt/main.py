# -*- coding: utf-8 -*-


""" Reads in crux/percolator tab-delimited results (psms) and returns tmt values"""

import os.path
import sys
import re
import pandas as pd
import tqdm
import logging

from pytmt import __version__
from pytmt.get_spec import Mzml
from pytmt import tmt_reporters

def quant(args):
    """
     reads in Percolator tab-delimited results (PSMS) \\
     and filter each row by protein-uniqueness and by Percolator q value \\
     then it opens the corresponding mzML file of the fraction and finds the scan \\
     and returns the intensity of each TMT reporter within specified \\
     mass tolerance. Currently supports only Percolator, and MS2-level quantification.

    Tested on:
        standalone comet 2017.01rev4 > crux 3.1 percolator
        crux 3.1 tide > crux 3.1 percolator

    To-do features:
        ms3 or multi-notch (shifting scan numbers)
        read in mzID files rather than percolator
        normalization and isotope purity adjustment
        filtering based on ms1

    Known issues:
        uses only directory index to match mzml files because of percolator

    Note:
        currently the sum of intensities is returned if multiple peaks are within the tolerance of reporter

    Usage:
        pytmt tests/data/mzml tests/data/percolator -o ~/Desktop/pytmt

    Example values for arguments:
        mzml_loc = 'tests/data/mzml'
        id_loc = 'tests/data/percolator'
        precision = 10
        q_filter = 0.1
        unique_only = True

    :param args:    arguments from argparse
    :return:        Exit OK
    """

    # Main logger setup
    main_log = logging.getLogger('pytmt')
    main_log.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    os.makedirs(args.out, exist_ok=True)
    fh = logging.FileHandler(os.path.join(args.out, 'tmt.log'))
    fh.setLevel(logging.INFO)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add the handlers to the logger
    main_log.addHandler(fh)
    main_log.addHandler(ch)

    main_log.info(args)
    main_log.info(__version__)

    # Folders of mzML files and the Percolator results.
    mzml_loc = args.mzml
    id_loc = args.id

    assert args.multiplex in [0, 2, 6, 10, 11, 16], '[error] TMT multiplexity not 0, 2, 6, 10, 11 or 16'

    # Get the reporter masses
    reporters = tmt_reporters.get_reporter(args.multiplex)

    # Define the PPM of integration
    precision = args.precision
    assert 1 <= precision <= 1000, '[error] mass tolerance must be between 1 and 1000 ppm'

    # Define the Percolator Q value cutoff filter
    q_filter = args.qvalue

    # Define the Percolator protein-unique PSM filter
    unique_only = args.unique

    assert os.path.isdir(mzml_loc), '[error] mzml directory not valid'
    assert os.path.isdir(id_loc), '[error] percolator directory not valid'

    # List all files in the percolator directory ending with target.psms.txt.
    id_files = [f for f in os.listdir(id_loc) if f.endswith('target.psms.txt')]

    assert len(id_files) == 1, 'Check percolator output directory has 1 *.target.psms.txt'

    # Read the Percolator psms file.
    id_df = pd.read_csv(filepath_or_buffer=os.path.join(id_loc, id_files[0]),
                        sep='\t')

    # Get all the file indices in the Percolator results file.
    file_indices = list(set(id_df['file_idx']))

    # Check that the number of mzMLs in the mzML folder is the same as the maximum of the ID file's file_idx column.
    # Note this will throw an error if not every fraction results in at least some ID, but we will ignore for now.

    mzml_files = [f for f in os.listdir(mzml_loc) if re.match('^.*.mzML', f)]

    # Sort the mzML files by names
    # Note this may create a problem if the OS Percolator runs on has natural sorting (xxxx_2 before xxxx_10)
    # But we will ignore for now
    mzml_files.sort()

    # Throw an error if there is no mzML file in the mzml directory
    assert len(mzml_files) != 0, '[error] no mzml files in the specified directory'
    assert len(mzml_files) == max(file_indices) + 1, '[error] number of mzml files not matching id list'

    # For each file index (fraction), open the mzML file, and create a subset Percolator ID dataframe
    for idx in file_indices:

        # Logging mzML
        main_log.info('Reading mzml file: {0} ({1} of {2})'.format(mzml_files[idx],
                                                                   str(idx + 1),
                                                                   str(len(file_indices)),
                                                                   )
                      )

        # Make a subset dataframe with the current file index (fraction) being considered
        fraction_id_df = id_df[id_df['file_idx'] == idx]

        # Arrange the PSM rows by scan number
        fraction_id_df = fraction_id_df.sort_values(by='scan').reset_index(drop=True)

        # Open the mzML file
        fraction_mzml = Mzml(path=os.path.join(mzml_loc, mzml_files[idx]),
                             precision=precision)
        fraction_mzml.parse_mzml_ms2()

        # Loop through each qualifying row in sub_df_filtered
        for i in tqdm.trange(len(fraction_id_df)):

            # Get current scan number
            scan = fraction_id_df.loc[i, 'scan']

            '''
            # Verbosity 1 progress message (display progress every 1000 rows)
            if (i+1) % 1000 == 0 and args.verbosity > 0:
                print('[verbosity 1] doing mzml:', mzml_files[idx], '(' + str(idx + 1),
                      'of', str(len(file_indices)) + ');', 'PSM:', i+1, 'of', len(fraction_id_df),
                      '(scan number:', str(scan) + ')', sep=' ')

            # Verbosity 2 progress message (calculates time remaining)
            if (i + 1) % 1000 == 0 and args.verbosity == 2:
                t2 = time()
                avg_time = (t2 - t1) / (i + 1)
                eta = ((len(fraction_id_df) - i) * avg_time) / 60

                print('[verbosity 2] elapsed:', round((t2 - t1) / 60, 2), 'minutes.',
                      'estimated remaining time for current fraction:', round(eta, 2), 'minutes.', sep=' ')
                      
            '''

            # If q-value filter is on, skip any row that fails the filter.
            if fraction_id_df.loc[i, 'percolator q-value'] > q_filter:
                continue

            # If the protein-unique PSM filter is on, skip any row that fails the filter.
            if unique_only and len(fraction_id_df.loc[i, 'protein id'].split(',')) > 1:
                continue

            # If this is a qualifying row, get the spectrum in the mzML file by scan number
            try:
                spectrum = fraction_mzml.msdata.get(scan)

            except KeyError:
                main_log.error('[error] spectrum index out of bound')
                continue

            #except xml.etree.ElementTree.ParseError:
            #    if args.verbosity == 2:
            #        print('[verbosity 2] XML eTree does not appear to be able to read this spectrum',
            #              '(scan number:', str(scan) + ')', sep=' ')
            #    continue

            #assert spectrum['ms level'] > 1, '[error] specified spectrum is a full scan'

            # For each reporter, check that the spectrum has that peak using pymzml
            # This returns a (m/z, intensity) tuple
            # We will append the intensity of each reporter to a list
            tmt_intensities = [idx, scan]

            for reporter in reporters:

                upper = reporter + reporter*(precision/2)*1e-6
                lower = reporter - reporter*(precision/2)*1e-6

                reporter_intensity = sum([I for mz_value, I in spectrum if upper > mz_value > lower])
                tmt_intensities.append(round(reporter_intensity, 2))

            # Write total spectrum intensity
            spectrum_intensity = sum([I for mz_value, I in spectrum])
            tmt_intensities.append(round(spectrum_intensity, 2))

            # Grow the results into an output list of lists (creating if does not exist)
            try:
                output_list.append(tmt_intensities)

            except NameError:
                output_list = [tmt_intensities]

    # Turn the output list into a data frame
    output_df_columns = ['file_idx', 'scan', ]

    for reporter in reporters:
        output_df_columns.append('m' + str(reporter))

    output_df_columns.append('spectrum_int')

    output_df = pd.DataFrame(output_list, columns=output_df_columns)

    # Final output, merging the input and output tables
    final_df = pd.merge(id_df, output_df, how='left')

    # Create output directory if it does not exist
    save_path = os.path.join(args.out, 'tmt_out.txt')

    # Save the file
    final_df.to_csv(save_path, sep='\t')

    main_log.info("Run completed successfully.")
    return sys.exit(os.EX_OK)


def main():
    """
    Entry point

    :return:
    """

    import argparse

    parser = argparse.ArgumentParser(description='pytmt returns ms2 tmt quantification'
                                                 'values from Crux Percolator output')

    parser.add_argument('mzml', help='path to folder containing mzml files')

    parser.add_argument('id', help='path to folder containing percolator tab-delimited files')

    parser.add_argument('-u', '--unique', action='store_true', help='quantify unique peptides only')

    parser.add_argument('-q', '--qvalue',
                        help='quantify peptides with q value below this threshold [default: 1.0]',
                        type=float,
                        default=1.0)

    parser.add_argument('-m', '--multiplex',
                        help='TMT-plex (0, 2, 6, 10, 11, 16) [default:10]',
                        type=int,
                        default=10)

    parser.add_argument('-p', '--precision',
                        help='ms2 spectrum mass shift tolerance in ppm [default: 10]',
                        type=int,
                        default=10)

    parser.add_argument('-o', '--out', help='name of the output directory [default: tmt_out]',
                        default='tmt_out')


    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    parser.set_defaults(func=quant)


    # Print help message if no arguments are given

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
