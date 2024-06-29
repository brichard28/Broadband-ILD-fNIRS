def raw_nirx_channelwise_fft(raw, plot, snr_thres):
    """

    :param raw:
    :param plot: bool
    :return:
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
    from matplotlib import pyplot as plt
    import matplotlib.ticker as ticker

    raw = raw.copy().load_data()
    freqs = np.unique(_channel_frequencies(raw.info))

    # only consider the long channels
    # raw_long = get_long_channels(raw)

    # picks = _check_channels_ordered(raw_long.info, freqs)
    picks = _check_channels_ordered(raw.info, freqs)

    # set up the figure with subplots for # channels * 2
    # n_channels = int(len(raw_long.ch_names) / 2)

    # need a list of channel locations
    # layout = find_layout(raw_long.info)
    layout = find_layout(raw.info)
    layout = deepcopy(layout)
    layout.pos[:, :2] -= layout.pos[:, :2].min(0)
    layout.pos[:, :2] /= layout.pos[:, :2].max(0)
    positions = layout.pos[:, :2] * 0.9

    # set up subplots
    fig = plt.figure(figsize=(5, 4), dpi=200)

    width, height = 0.05, 0.05

    SNR_dict = {}
    # snr_thres = 3.5

    # for each channel
    for ii in picks[::2]:

        pos = positions[ii, :]

        # take the fft of each channel, plot in a subplot - only look at 850 nm signal here
        # spectro_f, spectro_ch = signal.welch(raw_long._data[ii + 1, :], fs=raw.info['sfreq'])
        spectro_f, spectro_ch = signal.welch(raw._data[ii + 1, :], fs=raw.info['sfreq'])

        # convert to dB - referenced to minimum of each channel
        spectro_ch_db = 10*np.log10(spectro_ch / np.min(spectro_ch))

        # find the peak at the heart rate, then return the SNR value for each channel as a dictionary
        # looking between 45 and 96 bpm - this should capture > 99 % of all healthy resting people
        # - https://doi.org/10.1371/journal.pone.0227709.g002
        hr_peak_db = np.amax(spectro_ch_db[np.where((spectro_f > 0.8) & (spectro_f < 1.6))])
        hr_peak_db_ind = np.argmax(np.array(spectro_ch_db == hr_peak_db).astype(int))  # this might be bugged...

        # find the noise floor
        noise_floor = np.mean(spectro_ch_db[int(hr_peak_db_ind * 1.35):int(hr_peak_db_ind * 1.5)])

        # SNR_dict[f'{raw_long.ch_names[ii][:-4]}'] = hr_peak_db / spectro_ch_db[int(np.floor(hr_peak_db_ind * 1.5))]
        SNR_dict[f'{raw.ch_names[ii][:-4]}'] = hr_peak_db / noise_floor

        # if SNR_dict[f'{raw_long.ch_names[ii][:-4]}'] > snr_thres:
        if SNR_dict[f'{raw.ch_names[ii][:-4]}'] > snr_thres:
            cc = 'b'
        else:
            cc = 'r'

        # plot --- [lowerCorner_x, lowerCorner_y, width, height]
        ax = fig.add_axes([pos[0]+width/2, pos[1], width, height])
        ax.plot(spectro_f, spectro_ch_db, color=cc, linewidth=0.3)
        ax.set_ylim(bottom=0, top=35)
        ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=3, length=2, width=0.5, labelcolor='w')
        ax.patch.set_alpha(0.2)
        # ax.set_title(f'{raw_long.ch_names[ii][:-4]}', fontsize=3)
        ax.set_title(f'{raw.ch_names[ii][:-4]}', fontsize=3, pad=0)

    # add an empty plot with labels
    ax = fig.add_axes([0.05, 0.075, 1.5*width, 1.5*height])
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

    if plot:
        return fig, bad_channels_total, SNR_dict

    else:
        plt.close(fig)
        return 0, bad_channels_total, SNR_dict
