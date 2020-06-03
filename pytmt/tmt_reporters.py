# -*- coding: utf-8 -*-

""" Returns TMT reporter masses"""

all_reporters = [126.127726,  # 126
                 127.124761,  # 127N
                 127.131081,  # 127C
                 128.128116,  # 128N
                 128.134436,  # 128C
                 129.131471,  # 129N
                 129.137790,  # 129C
                 130.134825,  # 130N
                 130.141145,  # 130C
                 131.138180,  # 131N
                 131.144499,  # 131C (11-plex)
                 132.141535,  # 132N (Pro)
                 132.147855,  # 132C (Pro)
                 133.144890,  # 133N (Pro)
                 133.141210,  # 133C (Pro)
                 134.148245,  # 134C (Pro)
                 ]

def get_reporters(plex):
    """
    Define the reporter ion m/z values. Support Thermo TMT 0, 2, 6, 10, 11, 16-plex for now.
    :param plex: tandem mass tag multiplex type,
    :return: a list of accurate mass values to integrate
    :rtype: list
    """

    #

    if plex == 16:
        reporters = [i for n, i in enumerate(all_reporters)]

    elif plex == 11:
        reporters = [i for n, i in enumerate(all_reporters) if n in range(0, 11)]

    elif plex == 10:
        reporters = [i for n, i in enumerate(all_reporters) if n in range(0, 10)]

    elif plex == 6:
        reporters = [i for n, i in enumerate(all_reporters) if n in [0, 2, 4, 5, 8, 9]]

    elif plex == 2:
        reporters = [i for n, i in enumerate(all_reporters) if n in [0, 2]]

    elif plex == 0:
        reporters = [i for n, i in enumerate(all_reporters) if n in range(0, 1)]

    return reporters