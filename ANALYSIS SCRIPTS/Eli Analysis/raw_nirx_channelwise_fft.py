def raw_nirx_channelwise_fft(raw, plot, snr_thres=3):
    """
    Computes the channel-wise Signal-to-Noise Ratio (SNR) of heart rate signals using FFT.
    This function computes the SNR as the difference between the heart rate peak and the
    noise floor in dB. The heart rate peak is found between 0.75 and 1.6 Hz, and the
    noise floor is estimated around 1.5 times the heart rate peak frequency.

    :param raw: Raw MNE object containing the NIRS data.
    :param plot: bool, whether to plot the spectra (default is False).
    :param snr_thres: float, SNR threshold in dB
    :return: If `plot=True`, returns a figure, list of bad channels, and the SNR dictionary.
             If `plot=False`, returns `None`, list of bad channels, and the SNR dictionary.

    Author: Eli Bulger
    Date: March 2022
    """
    import numpy as np
    from mne.io import BaseRaw
    from mne.io.constants import FIFF
    from mne.utils import _validate_type, warn
    from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, \
        _check_channels_ordered
    from mne_nirs.channels import (get_long_channels,
                                   get_short_channels)
    from scipy import signal
    from mne.channels.layout import find_layout
    from copy import deepcopy
    from scipy.signal import welch
    from matplotlib import pyplot as plt
    import matplotlib.ticker as ticker

    raw = raw.copy().load_data()
    freqs = np.unique(_channel_frequencies(raw.info))
    picks = _check_channels_ordered(raw.info, freqs)

    layout = find_layout(raw.info)
    layout = deepcopy(layout)
    layout.pos[:, :2] -= layout.pos[:, :2].min(0)
    layout.pos[:, :2] /= layout.pos[:, :2].max(0)
    positions = layout.pos[:, :2] * 0.9

    if plot is True:
        fig = plt.figure(figsize=(5, 4), dpi=200)
        width, height = 0.05, 0.05

    SNR_dict = {}
    bad_channels = []

    # for each channel
    for ii in picks[::2]:
        pos = positions[ii, :]

        # Compute Welch spectrum
        spectro_f, spectro_ch = welch(raw._data[ii + 1, :], fs=raw.info['sfreq'])
        spectro_ch_db = 10 * np.log10(spectro_ch / np.mean(spectro_ch))

        # Identify heart rate peak (45-96 bpm = ~0.75-1.6 Hz)
        hr_indices = (spectro_f > 0.75) & (spectro_f < 1.6)
        hr_peak_db = np.max(spectro_ch_db[hr_indices])
        hr_peak_idx = np.argmax(spectro_ch_db[hr_indices])
        hr_peak_freq = spectro_f[hr_indices][hr_peak_idx]

        # Calculate noise floor around 1.5x heart rate frequency
        noise_center = 1.5 * hr_peak_freq
        noise_indices = (spectro_f > (noise_center - 0.05)) & (spectro_f < (noise_center + 0.05))
        noise_floor_db = np.mean(spectro_ch_db[noise_indices])

        # Compute SNR in dB
        SNR_db = hr_peak_db - noise_floor_db
        SNR_dict[raw.ch_names[ii][:-4]] = SNR_db

        snr_thresh_dB = 10*np.log10(snr_thres)

        if plot is True:
            # Determine if the channel passes the threshold
            if SNR_db < snr_thresh_dB:
                color = 'r'
            else:
                color = 'b'

            # Plot spectra
            ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])
            ax.plot(spectro_f, spectro_ch_db, color=color, linewidth=0.3)
            ax.set_ylim(bottom=0, top=35)
            ax.xaxis.set_major_locator(plt.MaxNLocator(2))
            ax.yaxis.set_major_locator(plt.MaxNLocator(2))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(labelsize=3, length=2, width=0.5, labelcolor='w')
            ax.patch.set_alpha(0.2)
            # ax.set_title(f'{raw_long.ch_names[ii][:-4]}', fontsize=3)
            ax.set_title(f'{raw.ch_names[ii][:-4]}', fontsize=3, pad=0)

    if plot is True:
        # add an empty plot with labels
        ax = fig.add_axes([0.05, 0.075, 1.5 * width, 1.5 * height])
        ax.plot(spectro_f, np.zeros((len(spectro_f),)), '-w', linewidth=0.3)
        ax.set_ylim(bottom=0, top=35)
        ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        ax.set_xlabel('Frequency (Hz)', fontsize=4)
        ax.set_ylabel('Power (dB)', fontsize=4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=4)

    # number of channels with relatively good SNR
    # n_good_channels = np.sum(np.array(list(SNR_dict.values())) > snr_thres)
    # good_chan_inds = [value > snr_thres for ch_name, value in SNR_dict.items()]
    # good_channels = np.array(list(SNR_dict.keys()))[good_chan_inds]

    bad_chan_inds = [value < snr_thres for ch_name, value in SNR_dict.items()]
    bad_channels = np.array(list(SNR_dict.keys()))[bad_chan_inds]

    bad_760 = [item + ' 760' for item in bad_channels]
    bad_850 = [item + ' 850' for item in bad_channels]

    bad_channels_total = bad_760 + bad_850

    if plot is True:
        return fig, bad_channels_total, SNR_dict

    else:
        return 0, bad_channels_total, SNR_dict
