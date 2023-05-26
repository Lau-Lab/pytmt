# -*- coding: utf-8 -*-


""" Perform matrix correction using a NNLS """

import io
import pandas as pd
import numpy as np
import scipy.optimize


def correct_matrix(output_df: pd.DataFrame,
                   contam: io.TextIOWrapper,
                   nnls: bool = True,
                   ) -> pd.DataFrame:
    """

    :param output_df:   TMT intensity output matrix, with file_idx,  scan, m..., spectrum_int as columns
    :param contam:      File handle to the contaminant matrix
    :param nnls:        Bool: uses non-negative least square for correction
    :return:            pd.Dataframe with Corrected TMT intensity output dataframe with additional columns

    """
    # Read the contaminant matrix
    contam = pd.read_csv(contam, index_col=0)

    # Normalize the contaminant matrix prior to correction
    contam_star = contam / contam.sum(axis=0)

    # Check dimension
    assert len(contam_star.columns) == len(output_df.columns) - 3, 'contaminant matrix not ' \
                                                                   'the same dimension as number of tags'

    if nnls:
        # Rowwise, use scipy NNLS to correct the output_df (TMT intensities) with the normalized contam matrix
        array_list = [scipy.optimize.nnls(contam_star, output_df.iloc[i, 2:len(output_df.columns) - 1])[0]
                  for i in range(0, len(output_df.index))]

    else:
        # Use numpy linalg solve instead
        array_list = [np.linalg.solve(contam_star, output_df.iloc[i, 2:len(output_df.columns) - 1])
                      for i in range(0, len(output_df.index))]

    new_column_names = [f'{rep}_cor' for rep in output_df.columns[2:len(output_df.columns)-1]]

    cor_df = pd.DataFrame(np.vstack(array_list), columns=new_column_names)

    # column bind the two data frames
    cor_output_df = pd.concat([output_df, cor_df], axis=1)

    return cor_output_df
