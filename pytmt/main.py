# -*- coding: utf-8 -*-

""" Reads in crux/percolator tab-delimited results (psms) and returns tmt values"""

import os.path
import sys
import re
import pandas as pd
import tqdm
import argparse

from pytmt import __version__
from pytmt.get_spec import Mzml
from pytmt import tmt_reporters
from pytmt import quantify_spec
from pytmt import correct_matrix

from pytmt.logger import get_logger
from pytmt.protein_group import get_canonical_parsimony_groups

def quant(args: argparse.Namespace) -> None:
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
        pytmt tests/data/mzml tests/data/percolator/percolator.target.psms.txt -o out

    Example values for arguments:
        mzml = 'tests/data/mzml'
        id = 'tests/data/percolator'
        precision = 10
        qvalue = 0.1

    :param args:    arguments from argparse
    :return:        Exit OK
    """

    # ---- Get the logger ----
    logger = get_logger(__name__, args.out)
    logger.info(args)
    logger.info(__version__)

    # Get the reporter masses
    reporters = tmt_reporters.get_reporters(args.multiplex)

    # Define the PPM of integration
    precision = args.precision

    assert os.path.isdir(args.mzml), '[error] mzml directory path not valid'
    assert os.path.isfile(args.id.name), '[error] percolator file path not valid'

    # List all files in the percolator directory ending with target.psms.txt.
        # id_files = [f for f in os.listdir(id_loc) if f.endswith('target.psms.txt')]
        # assert len(id_files) == 1, 'Check percolator output directory has 1 *.target.psms.txt'


    # Read the Percolator psms file.
    is_standalone = False
    try:
        id_df = pd.read_csv(filepath_or_buffer=args.id,  # os.path.join(id_loc, id_files[0]),
                            sep='\t')

    except pd.errors.EmptyDataError:
        logger.error('Unable to read percolator')
        sys.exit(1)

    # Check if the file is from standalone percolator or not

    except pd.errors.ParserError:
        is_standalone = True

    # If there is no file_idx column, then this is a standalone percolator file
    if is_standalone is False:
        try:
            id_df['file_idx']

        except KeyError:
            is_standalone = True

    if is_standalone:
        logger.info("Unable to find file_idx, attempting to read as standalone Percolator file")

        # The standalone percolator output has different number of columns per row because the proteins are separated by tabs
        # Read the first half of the table without the protein IDs
        with open(args.id.name, 'r') as f:
            f_ln = f.readlines()

        id_df = pd.DataFrame([ln.split('\t')[0:5] for ln in f_ln[1:]])
        id_df.columns = ['PSMId', 'score', 'percolator q-value', 'posterior_error_prob', 'peptide']
        id_df['percolator q-value'] = id_df['percolator q-value'].astype(float)
        id_df['posterior_error_prob'] = id_df['posterior_error_prob'].astype(float)

        # Create a sequence column for compatibility
        id_df['sequence'] = [pep[2:-2] for pep in id_df['peptide']]

        # Then read in the protein names and join them by comma instead of tab
        id_df['protein id'] = [','.join(ln.rstrip().split('\t')[5:]) for ln in f_ln[1:]]

        # Split the PSMId column to create file_idx, scan, and charge.
        # 2023-08-07 This now takes only the MSFragger output format which is in the format filename.scan.scan.charge_index
        id_df['charge'] = [psm.split('.')[-1] for psm in id_df['PSMId']]
        id_df['charge'] = [psm.split('_')[-2] for psm in id_df['charge']]
        id_df['charge'] = id_df['charge'].astype(int)
        id_df['scan'] = [psm.split('.')[-2] for psm in id_df['PSMId']]
        id_df['scan'] = id_df['scan'].astype(int)
        print(id_df['charge'])
        print(id_df['scan'])
        # The file name is the underscore('_') split until the last 3 parts, then rejoined by underscore
        # in case there are underscores in the filename. We then remove everything
        # We then remove all directories to get base name
        id_df['file_name'] = [os.path.basename('.'.join(psm.split('.')[:-3])) for psm in id_df['PSMId']]
        print(id_df['file_name'])

        # Get the sorted file names, hopefully this is the same index as the Crux Percolator output
        # TODO: Read the Percolator log file to get actual index and use file names to open the mzml instead
        sorted_index = sorted(set(id_df['file_name']))
        id_df['file_idx'] = id_df['file_name'].apply(sorted_index.index)

    # Get all the file indices in the Percolator results file.
    file_indices = list(set(id_df['file_idx']))

    # 2022-03-28 pytmt will now attempt to read the percolator.log.txt file for fraction (file_idx) mzML assignment
    log_path = os.path.join(os.path.dirname(args.id.name), 'percolator.log.txt')

    # If the log file exists, use it to read the assignment
    if os.path.exists(log_path):
        logger.warning(f'Percolator log file exists at {log_path} and will be used for index assignment.')

        with open(log_path, 'r') as f:
            lines = f.readlines()

        mzml_files = {}
        for line in lines:
            pattern = re.findall('INFO: Assigning index ([0-9]*) to (.*)\.', line)
            if len(pattern) == 1:
                idx, pathname = pattern[0]
                dirname, filename = os.path.split(pathname)
                mzml_files[int(idx)] = re.sub('\.pep\.xml', '', filename)
                # TODO: will probably have to account for .pin or other input to Percolator

    # If the log file does not exist, assign index naively based on sort
    else:
        logger.warning(f'Percolator log file not found at {log_path}; '
                         f'mzml files will be sorted for index assignment.')
        # Check that the number of mzMLs in the mzML folder is the same as the maximum of the ID file's file_idx column.
        # Note this will throw an error if not every fraction results in at least some ID, but we will ignore for now.

        mzml_filelist = [f for f in os.listdir(args.mzml) if re.match('^.*.mzML', f)]
        # Sort the mzML files by names
        # Note this may create a problem if the OS Percolator runs on has natural sorting (xxxx_2 before xxxx_10)
        # But we will ignore for now
        mzml_filelist.sort()

        # Make dictionary of idx, filename from enumerate
        mzml_files = {}
        for idx, filename in enumerate(mzml_filelist):
            mzml_files[idx] = re.sub('.mz[Mm][Ll](\.gz)?', '', filename)

        # Throw an error if there is no mzML file in the mzml directory
        assert len(mzml_files) != 0, '[error] no mzml files in the specified directory'
        assert len(mzml_files) == max(file_indices) + 1, '[error] number of mzml files not matching id list'

    # Print mzml files to log
    logger.info(f'mzml file orders: {mzml_files}')

    # For each file index (fraction), open the mzML file, and create a subset Percolator ID dataframe
    for idx in file_indices:

        # 2022-03-28 try to open either mzML or mzML.gz
        if os.path.exists(os.path.join(args.mzml, mzml_files[idx] + '.mzML')):
            mzml_path = os.path.join(args.mzml, mzml_files[idx] + '.mzML')
        elif os.path.exists(os.path.join(args.mzml, mzml_files[idx] + '.mzML.gz')):
            mzml_path = os.path.join(args.mzml, mzml_files[idx] + '.mzML.gz')
        else:
            raise FileNotFoundError(f'Could not find mzML file for index {idx} at {args.mzml}')

        # Logging mzML
        logger.info(f'Reading mzml file: {os.path.basename(mzml_path)} ({idx + 1} of {len(file_indices)})')

        # Make a subset dataframe with the current file index (fraction) being considered
        fraction_id_df = id_df[id_df['file_idx'] == idx]

        # Arrange the PSM rows by scan number
        fraction_id_df = fraction_id_df.sort_values(by='scan').reset_index(drop=True)

        # Open the mzML file
        fraction_mzml = Mzml(path=mzml_path,
                             precision=precision,
                             logger=logger,
                             )
        fraction_mzml.parse_mzml_ms2()
        if fraction_mzml.ms3data != {}:
            logger.info(f'Found {len(fraction_mzml.ms3data)} MS3 spectra in {os.path.basename(mzml_path)}')

        # Loop through each qualifying row in sub_df_filtered
        for i in tqdm.trange(len(fraction_id_df)):

            # Get current scan number
            scan = fraction_id_df.loc[i, 'scan']

            # If q-value filter is on, skip any row that fails the filter.
            if fraction_id_df.loc[i, 'percolator q-value'] > args.qvalue:
                continue

            # If the parsimony setting is set to unique, skip any row that fails the filter.
            if args.parsimony == 'unique' and len(fraction_id_df.loc[i, 'protein id'].split(',')) > 1:
                continue

            # If this is a qualifying row, get the spectrum in the mzML file by scan number
            # First check if the spectrum's ms3data is empty
            if fraction_mzml.ms3data == {}:  # if ms3data is empty
                try:
                    spectrum = fraction_mzml.ms2data.get(scan)
                    # Tidy this up
                    if spectrum is None:
                        logger.error(f'[error] spectrum index {scan} out of bound or is empty')
                        raise KeyError

                except KeyError:
                    logger.error('[error] spectrum index out of bound')
                    continue

            else:  # if ms3data is not empty
                try:
                    # Find the spectrum in the ms3data dictionary that has precursor id equal the ms2 spectrum
                    for key, value in fraction_mzml.prec_idx.items():
                        if int(value) == scan:
                            spectrum = fraction_mzml.ms3data.get(key)

                    # Tidy this up
                    if spectrum is None:
                        logger.error(f'[error] spectrum index {scan} out of bound or is empty')
                        raise KeyError

                except KeyError:
                    logger.error('[error] spectrum index out of bound')
                    continue




            # Get the intensity of each reporter
            tmt_intensities = quantify_spec.quantify_reporters(idx=idx,
                                                               scan=scan,
                                                               spectrum=spectrum,
                                                               precision=precision,
                                                               reporters=reporters,
                                                               digits=2,
                                                               )

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

    # Correct for contamination
    if args.contam is not None:
        output_df = correct_matrix.correct_matrix(output_df=output_df,
                                                  contam=args.contam,
                                                  nnls=args.nnls,
                                                  )

    # Final output, merging the input and output tables
    final_df = pd.merge(id_df, output_df, how='left')

    # Label light and heavy peptides
    if args.silac:
        heavy_mods = ["R\\[10.01\\]", "R\\[239.17\\]", "K\\[8.01\\]", "K\\[237.18\\]", "K\\[466.34\\]"]
        heavy = "|".join(heavy_mods)

        def add_heavy_tag(string: str) -> str:
            """ add _H to the end of each uniprot accession """
            return re.sub("(sp\\|)(.+)(\\|.*$)", "\\1\\2_H\\3", string)

        # Add _H to protein names if the peptide is heavy (contains the heavy tag)
        final_df['protein id'] = [",".join(map(add_heavy_tag, final_df['protein id'][i].split(',')))
                                  if bool(re.search(heavy, final_df['sequence'][i]))
                                  else final_df['protein id'][i]
                                  for i in range(len(final_df.index))]

    # Save the peptide file
    final_df.to_csv(os.path.join(args.out, 'tmt_out.txt'), sep='\t')

    # Collapse to protein level
    if args.contam is not None:
        protein_column_list = ['protein id', ] + [f'm{reporter}_cor' for reporter in reporters]
    else:
        protein_column_list = ['protein id'] + [f'm{reporter}' for reporter in reporters]

    protein_df = final_df[protein_column_list]

    # Sum the reporter intensities for each protein
    if args.parsimony == 'unique':
        filtered_protein_df = protein_df[protein_df['protein id'].str.count(',') == 0]

    elif args.parsimony == 'all':
        filtered_protein_df = protein_df

    elif args.parsimony == 'canonical':
        filtered_protein_df = get_canonical_parsimony_groups(result_df = protein_df,
                                                             contam=args.contam,
                                                             reporters=reporters,
                                                             )

    # Group by protein and sum the reporter intensities
    filtered_protein_df = filtered_protein_df.groupby('protein id').sum()

    # Remove any rows that are all zeros
    filtered_protein_df = filtered_protein_df[filtered_protein_df.sum(axis=1) > 0]

    # Save the protein file
    filtered_protein_df.to_csv(os.path.join(args.out, 'tmt_protein_out.txt'), sep='\t')

    logger.info("Run completed successfully.")

    return None #sys.exit(os.EX_OK)


class CheckReadableDir(argparse.Action):
    """ Class to check if directory is readable. """
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir = values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace, self.dest, prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


class CheckQValue(argparse.Action):
    """ Class to check that q values are between 0 and 1. """
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            values = float(values)
        except ValueError:
            raise argparse.ArgumentTypeError("%r for q_value not a floating-point literal" % (values,))

        if values < 0.0 or values > 1.0:
            raise argparse.ArgumentTypeError("%r for q_value not in range [0.0, 1.0]" % (values,))
        setattr(namespace, self.dest, values)


def main() -> None:
    """
    Entry point

    :return:
    """

    parser = argparse.ArgumentParser(description='pytmt returns ms2 tmt quantification values'
                                                 'from Percolator output and perform contamination'
                                                 'correction',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog='For more information, see GitHub repository at '
                                            'https://github.com/ed-lau/pytmt',
                                     )

    parser.add_argument('mzml',
                        help='<required> path to folder containing mzml files',
                        type=str,
                        action=CheckReadableDir,
                        )

    parser.add_argument('id',
                        help='<required> path to percolator target psms output file',
                        type=argparse.FileType('r'),
                        )

    parser.add_argument('-P', '--parsimony',
                        help='rule to collapse peptides into the protein level. '
                        'options: "all" (default; include all peptides), "unique"'
                             ' (only unique peptides), "canonical" (collapse'
                             'into parsimony groups).',
                        choices=['all', 'unique', 'canonical'],
                        default='all',
                        )

    parser.add_argument('-q', '--qvalue',
                        help='quantify peptides with q value below this threshold [default: 1.0]',
                        metavar='[0, 1]',
                        type=float,
                        default=1.0,
                        action=CheckQValue,
                        )

    parser.add_argument('-m', '--multiplex',
                        help='TMT-plex (0, 2, 6, 10, 11, 16, 18) [default:10]',
                        choices=[0, 2, 6, 10, 11, 16, 18],
                        type=int,
                        default=10)

    parser.add_argument('-p', '--precision',
                        help='ms2 spectrum mass shift tolerance in ppm [default: 10]',
                        type=int,
                        choices=[0, 100],
                        metavar='[0, 100]',
                        default=10)

    parser.add_argument('-o', '--out', help='name of the output directory [default: tmt_out]',
                        default='tmt_out')

    parser.add_argument('-c', '--contam',
                        help='Path to contaminant matrix csv file.'
                             ' Leave blank to get tmt output without correction',
                        type=argparse.FileType('r'),
                        )

    parser.add_argument('-n', '--nnls',
                        action='store_true',
                        help='uses non-negative least square for contamination correction',
                        )

    parser.add_argument('-S', '--silac',
                        help='mark peptides with SILAC reporter ions',
                        action='store_true',
                        )

    parser.add_argument('-v', '--version', action='version',
                        version='pyTMT {version}'.format(version=__version__))

    parser.set_defaults(func=quant)


    # Print help message if no arguments are given

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Convert args to a dictionary
    # args_dict = vars(args)

    # Run the function in the argument
    args.func(args)

    return None #sys.exit(os.EX_OK)


if __name__ == '__main__':
    main()