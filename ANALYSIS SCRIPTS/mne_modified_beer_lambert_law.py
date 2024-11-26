def mne_modified_beer_lambert_law(raw):
    """Convert NIRS optical density data to haemoglobin concentration.

    Parameters
    ----------
    raw : instance of Raw
        The optical density data.

    Returns
    -------
    raw : instance of Raw
        The modified raw instance.
    """

    import numpy as np
    from mne.io import BaseRaw
    from mne.io.constants import FIFF
    from mne.utils import _validate_type, warn
    from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, \
        _check_channels_ordered, _channel_chromophore, _validate_nirs_info

    raw = raw.copy().load_data()
    _validate_type(raw, BaseRaw, 'raw')
    # _validate_type(ppf, 'numeric', 'ppf')
    # ppf = float(ppf)
    # freqs = np.unique(_channel_frequencies(raw.info))
    # picks = _check_channels_ordered(raw.info, freqs)

    picks = _validate_nirs_info(raw.info, fnirs="od", which="Beer-lambert")
    freqs = np.array([raw.info["chs"][pick]["loc"][9] for pick in picks], float)

    abs_coef = _load_absorption(freqs)
    distances = source_detector_distances(raw.info)

    bad = ~np.isfinite(distances[picks])
    bad |= distances[picks] <= 0

    if bad.any():
        warn(
            "Source-detector distances are zero on NaN, some resulting "
            "concentrations will be zero. Consider setting a montage "
            "with raw.set_montage."
        )
    distances[picks[bad]] = 0.0
    if (distances[picks] > 0.1).any():
        warn(
            "Source-detector distances are greater than 10 cm. "
            "Large distances will result in invalid data, and are "
            "likely due to optode locations being stored in a "
            " unit other than meters."
        )
    rename = dict()

    # Modified beer lambert law equation
    """
    Solves the mBLL equation --> 
    __               __            __                                       __      __             __
    |    OD wl 1      |            |   DPF_1 . E_hbo_1    DPF_1 . E_hbr_1    |      |   d[C]_HbO    |
    |    OD wl 2      |    =       |   DPF_2 . E_hbo_2    DPF_2 . E_hbr_2    |  @   |   d[C]_HbR    |
    __               __            __                                       __      __             __
    """

    # find the DPF
    dpf760 = 6.2966  # from NIRx aurora manual
    dpf850 = 5.23433
    # coef_oxy_760 = 1.3495580e-01;       # extinction coefficient   mm - 1 / mM
    # coef_oxy_850 = 2.4365740e-01;
    # coef_deoxy_760 = 3.5662416e-01;
    # coef_deoxy_850 = 1.5921100e-01;

    """
     forwardmat = [dpf760*coef_oxy_760*dist(i), dpf760*coef_deoxy_760*dist(i); 
                    dpf850*coef_oxy_850*dist(i), dpf850*coef_deoxy_850*dist(i)];
    hemoglobin = 1000*[dOD(:, ch760(i)), dOD(:, ch850(i))] * pinv(forwardmat)'; % uM
    HBO(:,i) = hemoglobin(:,1);
    HBD(:,i) = hemoglobin(:,2);
    """

    picks_760 = np.arange(len(raw.ch_names)/2, dtype=int)
    picks_850 = np.arange(len(raw.ch_names)/2, len(raw.ch_names), dtype=int)

    # code the distances too

    # for ii in picks[::2]:
    # for ii, jj in zip(picks[::2], picks[1::2]):
    for ii, jj in zip(picks_760, picks_850):

        # wavelengths must be 760 and 850
        forwardmat = np.reshape(np.array([abs_coef[0, 0] * dpf760 * distances[ii], abs_coef[0, 1]* dpf760 * distances[ii],
                                  abs_coef[1, 0] * dpf850 * distances[ii], abs_coef[1, 1] * dpf850 * distances[ii]]), (2, 2))

        inverse_mat = np.linalg.inv(forwardmat)

        raw._data[[ii, jj]] = inverse_mat @ raw._data[[ii, jj]] * 1e-3

        # Update channel information
        coil_dict = dict(hbo=FIFF.FIFFV_COIL_FNIRS_HBO, hbr=FIFF.FIFFV_COIL_FNIRS_HBR)
        for ki, kind in zip((ii, jj), ("hbo", "hbr")):
            ch = raw.info["chs"][ki]
            ch.update(coil_type=coil_dict[kind], unit=FIFF.FIFF_UNIT_MOL)
            new_name = f'{ch["ch_name"].split(" ")[0]} {kind}'
            rename[ch["ch_name"]] = new_name
    raw.rename_channels(rename)

    # Validate the format of data after transformation is valid
    _validate_nirs_info(raw.info, fnirs="hb")

    return raw


def _load_absorption(freqs):
    """Load molar extinction coefficients."""
    # Data from https://omlc.org/spectra/hemoglobin/summary.html
    # The text was copied to a text file. The text before and
    # after the table was deleted. The the following was run in
    # matlab
    # extinct_coef=importdata('extinction_coef.txt')
    # save('extinction_coef.mat', 'extinct_coef')
    #
    # Returns data as [[HbO2(freq1), Hb(freq1)],
    #                  [HbO2(freq2), Hb(freq2)]]

    import os.path as op
    import numpy as np
    from scipy.io import loadmat
    from scipy.interpolate import interp1d

    extinction_fname = op.join(r'C:\Users\benri\Anaconda3\envs', 'spyder-env', 'Lib', 'site-packages',
                               'mne', 'data',
                               'extinction_coef.mat')
    a = loadmat(extinction_fname)['extinct_coef']

    interp_hbo = interp1d(a[:, 0], a[:, 1], kind='linear')
    interp_hb = interp1d(a[:, 0], a[:, 2], kind='linear')

    ext_coef = np.array([[interp_hbo(freqs[0]), interp_hb(freqs[0])],
                         [interp_hbo(freqs[1]), interp_hb(freqs[1])]])
    abs_coef = ext_coef * 0.2303

    return abs_coef
