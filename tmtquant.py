"""
TMT quantifier v.0.1.0
Reads in Crux/Percolator tab-delimited results (PSMS) and returns TMT values

Currently supports only Percolator, and MS2 quant.
To do: TMT 6-plex, MS3 or multi-notch, etc.

Molecular Proteomics Laboratory
http://maggielab.org

"""

import pymzml
import os.path
import re
import pandas as pd
import xml

# Folders of mzML files and the Percolator results.
mzml_loc = './test_mzml_2'
id_loc = './test_perc_2'

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
precision = 20

# Define the Percolator Q value cutoff filter
qFilter = 0.05

# Define the Percolator protein-unique PSM filter
uniqueOnly = True

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

    # Loop through each qualifying row in sub_df_filtered
    for i in range(len(fraction_id_df)):

        # If this is a qualifying PSM, pull the spectrum based on the scan number
        scan = fraction_id_df.loc[i, 'scan']

        if (i+1) % 100 == 0:
            print('Doing mzML:', mzml_files[idx], '(' + str(idx + 1), 'of', str(len(file_indices)) + ');',
                  'PSM:', i+1, 'of', len(fraction_id_df), '(scan number:', str(scan) + ')', sep=' ')

        # If q-value filter is on, skip any row that fails the filter.
        if fraction_id_df.loc[i, 'percolator q-value'] > qFilter:
            continue

        # If the protein-unique PSM filter is on, skip any row that fails the filter.
        if uniqueOnly and len(fraction_id_df.loc[i, 'protein id'].split(',')) > 1:
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



final_df.to_csv('output.csv')