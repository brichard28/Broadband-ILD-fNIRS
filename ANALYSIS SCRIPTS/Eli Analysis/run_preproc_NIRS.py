"""

Eli Bulger

March 2022


"""


def preprocess_NIRX(data, data_snirf=0, event_dict=0,
                    save=False, savename=None,
                    plot_steps=False,
                    crop=True, crop_low=None, crop_high=None,
                    events_from_snirf=True, events_modification=None,
                    drop_short=False, reject=True,
                    short_regression=True,
                    negative_enhancement=False,
                    snr_thres=3.5,
                    sci_thres=0.8,
                    filter_type='iir', filter_limits=[0.01, 0.3], # 0.01, 0.1
                    tddr=True):
    """

    :param data:
    :param data_snirf:
    :param event_dict:
    :param save:
    :param savename:
    :param plot_steps:
    :param crop:
    :param crop_low:
    :param crop_high:
    :param events_from_snirf:
    :param events_modification:
    :param drop_short:
    :param reject:
    :param short_regression:
    :param negative_enhancement:
    :param snr_thres:
    :param filter_type:
    :param filter_limits:
    :param tddr:
    :return:
    """

    import mne
    import mne_nirs
    from matplotlib import pyplot as plt
    import itertools
    from mne.preprocessing.nirs import (temporal_derivative_distribution_repair)
    from mne_modified_beer_lambert_law import mne_modified_beer_lambert_law
    #from mne_short_channel_correction import short_channel_regression_OD
    from raw_nirx_channelwise_fft import raw_nirx_channelwise_fft
    import numpy as np

    # ---------------------------------------------------------------
    # -----------------      Crop The Data        ---------
    # ---------------------------------------------------------------
    if crop is True:
        data_snirf_crop = data_snirf.crop(tmin=crop_low, tmax=crop_high)
        data_crop = data.crop(tmin=crop_low, tmax=crop_high)
    else:
        data_crop = data
        data_snirf_crop = data_snirf

    # ---------------------------------------------------------------
    # -----------------      Find Events        ---------
    # ---------------------------------------------------------------
    if events_from_snirf:
        events, event_dict2 = mne.events_from_annotations(data_snirf_crop, event_dict)
    else:
        events, event_dict2 = mne.events_from_annotations(data_crop, event_dict)

    # if events are doubled
    # if events_modification:
    #     events = events[1:-1:2, :]

    trial_annotations = mne.annotations_from_events(events, sfreq=data_snirf.info['sfreq'])
    data.set_annotations(trial_annotations)

    # ---------------------------------------------------------------
    # -----------------   Plot                  ---------
    # ---------------------------------------------------------------
    if plot_steps:
        x = 1
        # data_snirf_crop.plot(n_channels=30,
        #                      show_scrollbars=True,
        #                      scalings='auto')

    # ---------------------------------------------------------------
    # -----------------   conversion to OD                  ---------
    # ---------------------------------------------------------------
    raw_od = mne.preprocessing.nirs.optical_density(data_snirf_crop)

    if plot_steps:
        x = 1
        # raw_od.plot(n_channels=30,
        #             show_scrollbars=True,
        #             scalings='auto')

    # ---------------------------------------------------------------
    # -----------------   Scalp Coupling Index              ---------
    # ---------------------------------------------------------------
    if reject is True:
        sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)
        #fig, ax = plt.subplots()
        #ax.hist(sci)
        #ax.set(xlabel='Scalp Coupling Index', ylabel='Count', xlim=[0, 1])

        # raw_od.info['bads'] = list(itertools.compress(raw_od.ch_names, sci < 0.8))

        # visually inspect and mark bads at this point...
        # raw_od.plot(show_scrollbars=True,
        #             scalings='auto')

        # ---------------------------------------------------------------
        # -----------------   Plot the FFT of each channel      ---------
        # ---------------------------------------------------------------
        # rejecting based on SNR from the fft plot - an alternative method
        fig_fft, bad_channels_total, SNR_dict = raw_nirx_channelwise_fft(raw_od, plot_steps, snr_thres=snr_thres)
        # raw_od.info['bads'] = bad_channels_total

        agree = np.intersect1d(list(itertools.compress(raw_od.ch_names, sci < sci_thres)), bad_channels_total)
        diff = np.setdiff1d(list(itertools.compress(raw_od.ch_names, sci < sci_thres)), bad_channels_total)

        # set the bad channels to what the two methods agree upon
        raw_od.info['bads'] = list(np.concatenate([agree, diff]))
        #raw_od.info['bads'] = list(itertools.compress(raw_od.ch_names, sci < sci_thres)) + list(bad_channels_total)


    # ---------------------------------------------------------------
    # ----      Temporal derivative filter        -------
    # ---------------------------------------------------------------
    if tddr is True:
        corrected_tddr = temporal_derivative_distribution_repair(raw_od)
    else:
        corrected_tddr = raw_od

    # if plot_steps:
    #     corrected_tddr.plot(n_channels=30,
    #                         show_scrollbars=True,
    #                         scalings='auto')

    # ---------------------------------------------------------------
    # -----               High-pass Filter the Data            ------
    # ---------------------------------------------------------------
    if filter_limits[0] is not None:
        if filter_type == 'iir':
            raw_od_corrected = corrected_tddr.filter(l_freq=filter_limits[0], h_freq=None,
                                                      l_trans_bandwidth=filter_limits[0],
                                                      method='iir', phase='zero',
                                                      iir_params={'order': 3, 'ftype': 'butter', 'output': 'sos'})

        elif filter_type == 'fir':
            raw_od_corrected = corrected_tddr.filter(l_freq=filter_limits[0], h_freq=None,
                                                      l_trans_bandwidth=filter_limits[0],
                                                      method='fir')
    else:
        raw_od_corrected = corrected_tddr

    raw_od_corrected = corrected_tddr
    # ---------------------------------------------------------------
    # ----      Remove the short-channel information          -------
    # ---------------------------------------------------------------
    if short_regression is True:
        raw_OD_sc = mne_nirs.signal_enhancement.short_channel_regression(raw_od_corrected)
        # raw_OD_sc = short_channel_regression_OD(raw_od_corrected)
    else:
        raw_OD_sc = raw_od_corrected

    # if plot_steps:
    #     filter_fig, filter_axes = plt.subplots(nrows=2,ncols=1)
    #     curr_ax = filter_axes[0]
    #     raw_OD_sc.plot_psd(average='mean', ax = filter_axes[0])
    #     curr_ax.set_xlim([0, 0.3])
    #     curr_ax.set_title('Before filtering', weight='bold', size='x-large')

    # # ---------------------------------------------------------------
    # # -----               Low-pass Filter the Data            ------
    # # ---------------------------------------------------------------
    if filter_limits[1] is not None:
        if filter_type == 'iir':
            raw_OD_filt = raw_OD_sc.filter(l_freq=None, h_freq=filter_limits[1],
                                                 h_trans_bandwidth=filter_limits[1] / 2, method='iir', phase='zero',
                                                 iir_params={'order': 3, 'ftype': 'butter', 'output': 'sos'})

        elif filter_type == 'fir':
            raw_OD_filt = raw_OD_sc.filter(l_freq=None, h_freq=filter_limits[1],
                                                 h_trans_bandwidth=filter_limits[1] / 2, method='fir', phase='zero')
    else:
        raw_OD_filt = raw_OD_sc

    # if plot_steps:
    #     curr_ax = filter_axes[1]
    #     raw_OD_filt.plot_psd(average='mean',ax = curr_ax)
    #     curr_ax.set_xlim([0, 0.3])
    #     curr_ax.set_title('After filtering', weight='bold', size='x-large')
        

    # if plot_steps:
    #     raw_OD_filt.plot(n_channels=30,
    #                         show_scrollbars=True,
    #                         scalings='auto')

    # ---------------------------------------------------------------
    # ----      conversion to Hb and remove heart rate        -------
    # ---------------------------------------------------------------
    raw_haemo_filt = mne_modified_beer_lambert_law(raw_OD_filt)

    # if plot_steps:
    #     raw_haemo_filt.plot(n_channels=30,
    #                       show_scrollbars=True,
    #                       scalings='auto')
    #     fig = raw_haemo_filt.plot_psd(average=True)
    #     fig.suptitle('Before filtering', weight='bold', size='x-large')
    #     fig.subplots_adjust(top=0.88)

    # ---------------------------------------------------------------
    # -----              Plot and Drop the Short Channels            ------
    # ---------------------------------------------------------------
    if drop_short is True:
        # get only the long channels now
        raw_haemo_filt = mne_nirs.channels.get_long_channels(raw_haemo_filt)


    # ---------------------------------------------------------------
    # ------            Negative Enhancement                ---------
    # ---------------------------------------------------------------
    # Makes HbO and HbR more negatively correlated
    if negative_enhancement is True:
        raw_haemo_filt = mne_nirs.signal_enhancement.enhance_negative_correlation(raw_haemo_filt)

    # ---------------------------------------------------------------
    # -----------------          Save         ---------
    # ---------------------------------------------------------------
    # fix the event times bc MNE is bugged out >:(
    # events[:, 0] = events[:, 0] - crop_low * raw_haemo_filt.info['sfreq']

    if save is True:
        raw_haemo_filt.save(savename, overwrite=True)

    return raw_haemo_filt, events
