def mark_aux_movement_bad(raw, aux_snirf, cooldown_dur=5):
    """

    :param raw:
    :param aux_snirf:
    :param cooldown_dur:
    :return:
    """

    import numpy as np
    from scipy import stats

    aux_times = np.array(aux_snirf.axes[0].T)
    aux_fs = 1 / (aux_times[3] - aux_times[2])

    # convert all the metrics into zscores
    a1 = np.array(aux_snirf.accelerometer_1_x)
    a2 = np.array(aux_snirf.accelerometer_1_y)
    a3 = np.array(aux_snirf.accelerometer_1_z)
    g1 = np.array(aux_snirf.gyroscope_1_x)
    g2 = np.array(aux_snirf.gyroscope_1_y)
    g3 = np.array(aux_snirf.gyroscope_1_z)

    movement_metrics = np.stack([a1, a2, a3, g1, g2, g3])
    movement_avg = stats.zscore(np.mean(movement_metrics, axis=0))

    # find the locations where there is excess movement - greater than 2 or 3 standard deviations from the mean
    movement_excess = movement_avg > 3

    aux_onset = []

    # lame way but /e
    for i in range(len(movement_excess) - 1):
        if movement_excess[i] == False and movement_excess[i + 1] == True:
            # mark the index of all the False to True sections
            aux_onset.append(aux_times[i])

    # set the duration after each onset to be fairly long
    # (say, 5 seconds - based on qualitative analysis of time to reach baseline)
    aux_durations = np.ones((len(aux_onset),)) * cooldown_dur

    raw_annot = raw.copy()

    raw_annot.annotations.append(onset=aux_onset,
                           duration=aux_durations,
                           description=np.repeat(['BAD_movement'], len(aux_durations)),
                           ch_names=None)  # follow mne convention to reject

    return raw_annot
