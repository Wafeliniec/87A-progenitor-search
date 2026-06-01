"""
Gaia DR3 filtering function adapted from the data-reproduction
package accompanying Stoop et al. (2024), Nature, "Two waves of
massive stars running away from the young cluster R136".

Used here for a Gaia DR3 search around SN1987A.
Modified by: Lena Bukowska
Modifications: file paths/project context only; filtering logic initially unchanged.
"""
import numpy as np
import math
import pandas as ps
import matplotlib.pyplot as plt
from zero_point import zpt

def gaia_filtering(data, max_ruwe, max_ipd_frac_multi_peak, min_visibility_periods_used, max_duplicated_source, max_ipd_gof_harmonic_amplitude):
    """
    Applies corrections and filters to the raw Gaia DR3 data
    Input:
	data as pandas dataframe
	cut-off values for gaia parameters
    Output:
	corrected and filtered data as pandas dataframe
    """

    # correct for parallax zero-point offset
    zpt.load_tables()
    zero_point_5parm = zpt.get_zpt(data.phot_g_mean_mag.to_numpy()[data.astrometric_params_solved == 31], data.nu_eff_used_in_astrometry.to_numpy()[data.astrometric_params_solved == 31], data.pseudocolour.to_numpy()[data.astrometric_params_solved == 31], data.ecl_lat.to_numpy()[data.astrometric_params_solved == 31], data.astrometric_params_solved.to_numpy()[data.astrometric_params_solved == 31])
    zero_point_6parm = zpt.get_zpt(data.phot_g_mean_mag.to_numpy()[data.astrometric_params_solved == 95], data.nu_eff_used_in_astrometry.to_numpy()[data.astrometric_params_solved == 95], data.pseudocolour.to_numpy()[data.astrometric_params_solved == 95], data.ecl_lat.to_numpy()[data.astrometric_params_solved == 95], data.astrometric_params_solved.to_numpy()[data.astrometric_params_solved == 95])

    # zero-point offset correction needs to be SUBSTRACTED
    data.loc[data['astrometric_params_solved'] == 31, 'parallax'] -= zero_point_5parm
    data.loc[data['astrometric_params_solved'] == 95, 'parallax'] -= zero_point_6parm

    # remove sources with too high ruwe
    data = data[(data.ruwe < max_ruwe)|(data.astrometric_params_solved == 3)]
    print('After excluding on ruwe: %i sources left'%len(data))

    # remove sources with too high ipd_frac_multi_peak
    data = data[data.ipd_frac_multi_peak <= max_ipd_frac_multi_peak]
    print('After excluding on ipd_frac_multi_peak: %i sources left'%len(data))

    # remove sources with too low visibility_periods_used
    data = data[data.visibility_periods_used >= min_visibility_periods_used]
    print('After excluding on visibility_periods_used: %i sources left'%len(data))

    # remove sources with duplicated source
    data = data[data.duplicated_source < max_duplicated_source]
    print('After excluding on duplicated_source: %i sources left'%len(data))

    # remove sources with too high ipd_gof_harmonic_amplitude
    data = data[data.ipd_gof_harmonic_amplitude < max_ipd_gof_harmonic_amplitude]
    print('After excluding on ipd_gof_harmonic_amplitude: %i sources left'%len(data))

    return data
