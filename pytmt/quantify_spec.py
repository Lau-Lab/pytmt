# -*- coding: utf-8 -*-


""" Given a spectrum, precision, and list of reporters, get reporter intensity values """


def quantify_reporters(idx, scan, spectrum, precision, reporters, digits):
    """

    :param idx: file index, for reporting only
    :param scan: scan number of the spectrum, for reporting only
    :param spectrum: spectrum list of [mz/I]
    :param precision: mass precision
    :param reporters: list of reporters to be quantified
    :param digits: number of significant digits to report
    :return: list of intensities

    """

    tmt_intensities = [idx, scan]

    for reporter in reporters:
        upper = reporter + reporter * (precision / 2) * 1e-6
        lower = reporter - reporter * (precision / 2) * 1e-6

        try:
            reporter_intensity = sum([I for mz_value, I in spectrum if upper > mz_value > lower])
        except TypeError:
            reporter_intensity = 0

        tmt_intensities.append(round(reporter_intensity, digits))

    # Write total spectrum intensity
    spectrum_intensity = sum([I for mz_value, I in spectrum])
    tmt_intensities.append(round(spectrum_intensity, digits))

    return tmt_intensities
