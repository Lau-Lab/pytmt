"""
TMT-quant v.0.1.0
Reads in Crux/Percolator tab-delimited results (PSMS) and returns TMT values


Molecular Proteomics Laboratory
http://maggielab.org

"""

import pymzml
import os.path
import re
import pandas as pd
from time import time
import xml

def quant(args):
    """
    Reads in Crux/Percolator tab-delimited results (PSMS), open mzML files, and returns TMT values

    Currently supports only Percolator, and MS2 quant.
    To-do features: TMT 6-plex, MS3 or multi-notch, read in mzID, etc.

    Notes:
    mzml_loc = './test_mzml_2'
    id_loc = './test_perc_2'
    precision = 20
    qFilter = 0.05
    uniqueOnly = True

    :param args:    arguments from argaprse
    :return:
    """

    # Folders of mzML files and the Percolator results.
    mzml_loc = args.mzml
    id_loc = args.id

    # Define the reporter ion m/z values. Support TMT-10plex for now.
    reporters = [126.127726,
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

    # Define the PPM of integration
    precision = args.precision

    # Define the Percolator Q value cutoff filter
    qfilter = args.qvalue

    # Define the Percolator protein-unique PSM filter
    uniqueonly = args.unique

    assert os.path.isdir(mzml_loc), 'Percolator directory not valid.'
    assert os.path.isdir(id_loc), 'Percolator directory not valid.'

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
    assert len(mzml_files) != 0, 'No mzML files in the specified directory.'
    assert len(mzml_files) == max(file_indices) + 1, 'Number of files not matching.'



    # For each file index (fraction), open the mzML file, and create a subset Percolator ID dataframe
    for idx in file_indices:

        print('Doing mzML:', mzml_files[idx], '(' + str(idx + 1), 'of', str(len(file_indices)) + ')', sep= ' ')

        # Make a subset dataframe with the current file index (fraction) being considered
        fraction_id_df = id_df[id_df['file_idx'] == idx]

        # Arrange the PSM rows by scan number
        fraction_id_df = fraction_id_df.sort_values(by='scan').reset_index(drop=True)

        # Open the mzML file
        fraction_mzml = pymzml.run.Reader(path=os.path.join(mzml_loc, mzml_files[idx]),
                                      MS1_Precision=precision*1e-6, MSn_Precision=precision*1e-6)

        # Time counter at the start of the run
        t1 = time()

        # Loop through each qualifying row in sub_df_filtered
        for i in range(len(fraction_id_df)):

            # If this is a qualifying PSM, pull the spectrum based on the scan number
            scan = fraction_id_df.loc[i, 'scan']

            # Verbosity 1 progress message
            if (i+1) % 100 == 0 and args.verbosity > 0:
                print('verbosity 1: doing mzML:', mzml_files[idx], '(' + str(idx + 1),
                      'of', str(len(file_indices)) + ');', 'PSM:', i+1, 'of', len(fraction_id_df),
                      '(scan number:', str(scan) + ')', sep=' ')

            # Verbosity 2 progress message
            if (i + 1) % 100 == 0 and args.verbosity == 2:
                t2 = time()
                avg_time = (t2 - t1) / (i + 1)
                eta = ((len(fraction_id_df) - i) * avg_time) / 60

                print('verbosity 2: elapsed:', round((t2 - t1) / 60, 2), 'minutes.',
                      'estimated remaining time for fraction:', eta, 'minutes.', sep=' ')

            # If q-value filter is on, skip any row that fails the filter.
            if fraction_id_df.loc[i, 'percolator q-value'] > qfilter:
                continue

            # If the protein-unique PSM filter is on, skip any row that fails the filter.
            if uniqueonly and len(fraction_id_df.loc[i, 'protein id'].split(',')) > 1:
                continue

            # Get the spectrum in the mzML file by scan number
            try:
                spectrum = fraction_mzml[scan]

            except KeyError:
                print('Spectrum index out of bound.')

            except xml.etree.ElementTree.ParseError:
                print('XML eTree does not appear to be able to read this spectrum')
                continue

            assert spectrum['ms level'] > 1, 'Spectrum is not MSn.'

            # For each reporter, check that the spectrum has that peak using pymzml
            # This returns a m/z and intensity tuple
            # We will append the intensity of each reporter to a list
            tmt_intensities = [idx, scan]

            for reporter in reporters:
                matchList = spectrum.has_peak(reporter)

                if matchList and len(matchList) == 1:
                    for mz, I in matchList:
                        tmt_intensities.append(I)
                else:
                    tmt_intensities.append(0)

            # Grow the results into an output list of lists (creating if does not exist)
            try:
                output_list.append(tmt_intensities)
            except NameError:
                output_list = [tmt_intensities]

    # Turn the output list into a data frame
    output_df_columns = ['file_idx', 'scan', ]

    for reporter in reporters:
        output_df_columns.append('m' + str(reporter))

    output_df = pd.DataFrame(output_list, columns=output_df_columns)

    # Final output, merging the input and output tables
    final_df = pd.merge(id_df, output_df, how='left')

    # Create output directory if it does not exist
    os.makedirs(args.out, exist_ok=True)
    save_path = os.path.join(args.out, 'tmt_out.txt')

    final_df.to_csv(save_path, sep='\t'))

    return sys.exit(os.EX_OK)


"""
Argparse code for running main function with parsed arguments from command line
"""
if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='TMT-quant returns TMT quantification'
                                                 'values from Percolator output')

    parser.add_argument('mzml', help='path to folder containing mzml files')

    parser.add_argument('id', help='path to folder containing percolator tab-delimited files')

    parser.add_argument('-u', '--unique', action='store_true', help='quantify unique peptides only')

    parser.add_argument('-q', '--qvalue',
                        help='quantify peptides with q value below this threshold [default: 1]',
                        type=float,
                        default=1.0)

    parser.add_argument('-p', '--precision',
                        help='ms2 spectrum mass shift tolerance in ppm [default: 20]',
                        type=int,
                        default=20)

    parser.add_argument('-o', '--out', help='name of the output directory [default: tmt_out]',
                        default='tmt_out')

    parser.add_argument('-v', '--verbosity',
                        help='verbosity of error messages. 0=quiet, 1=default, 2=verbose [default: 1]',
                        type=int,
                        choices=[0, 1, 2],
                        default=1)

    parser.set_defaults(func=quant)

    # Print help message if no arguments are given
    import sys
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)

