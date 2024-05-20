



# def plot_nirs_evoked(fig, evoked, colour, nirs_pick):
#     """
#
#     :param fig:
#     :param evoked:
#     :param colour:
#     :param nirs_pick:
#     :return:
#     """
#
#     import matplotlib
#     from matplotlib import pyplot as plt
#     from mne.channels.layout import find_layout
#     from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, \
#         _check_channels_ordered
#     import numpy as np
#     from copy import deepcopy
#
#     layout = find_layout(evoked.info)
#     layout = deepcopy(layout)
#     layout.pos[:, :2] -= layout.pos[:, :2].min(0)
#     layout.pos[:, :2] /= layout.pos[:, :2].max(0)
#     positions = layout.pos[:, :2] * 0.9
#
#     freqs = np.unique(_channel_frequencies(evoked.info, nominal=True))
#     picks = _check_channels_ordered(evoked.info, freqs)
#
#     width, height = 0.05, 0.05
#
#     for ii in picks[::2]:
#
#         pos = positions[ii, :]
#
#         # plot --- [lowerCorner_x, lowerCorner_y, width, height]
#         ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])
#         ax.plot(evoked.times(), evoked.data[ii, :], color=colour, linewidth=0.3)
#         ax.set_ylim(bottom=0, top=35)
#         ax.xaxis.set_major_locator(plt.MaxNLocator(2))
#         ax.yaxis.set_major_locator(plt.MaxNLocator(2))
#         ax.spines['top'].set_visible(False)
#         ax.spines['right'].set_visible(False)
#         ax.tick_params(labelsize=3, length=2, width=0.5, labelcolor='w')
#         ax.patch.set_alpha(0.2)
#         # ax.set_title(f'{raw_long.ch_names[ii][:-4]}', fontsize=3)
#         ax.set_title(f'{evoked.ch_names[ii][:-4]}', fontsize=3, pad=0)
#
#     # add a legend
#     if nirs_pick == 'ΔHbO':
#         leg_lines = [line for line in ax.lines if line.get_c() == 'r'][:1]
#     else:
#         leg_lines = [line for line in ax.lines if line.get_c() == 'b'][:1]
#     fig.legend(leg_lines, f'{nirs_pick}', loc='lower right')
#
#     # add an empty plot with labels
#     ax = fig.add_axes([0.05, 0.075, 1.5 * width, 1.5 * height])
#     ax.plot(evoked.times(), np.zeros((len(evoked.times()),)), '-w', linewidth=0.3)
#     ax.set_ylim(bottom=0, top=35)
#     ax.xaxis.set_major_locator(plt.MaxNLocator(2))
#     ax.yaxis.set_major_locator(plt.MaxNLocator(2))
#     ax.set_xlabel('Frequency (Hz)', fontsize=4)
#     ax.set_ylabel('Power (dB)', fontsize=4)
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.tick_params(labelsize=4)
#
#     return fig



def plot_nirs_evoked_error(fig, evoked, error, colour, ylim, add_legend=False, hbdiff=False):
    """

    :param fig:
    :param evoked:
    :param error:
    :param colour:
    :param ylim:
    :param add_legend:
    :return:
    """

    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from mne.channels.layout import find_layout
    from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, \
        _check_channels_ordered
    import numpy as np
    from copy import deepcopy
    from mne.io.pick import _picks_to_idx

    layout = find_layout(evoked.info)
    layout = deepcopy(layout)
    layout.pos[:, :2] -= layout.pos[:, :2].min(0)
    layout.pos[:, :2] /= layout.pos[:, :2].max(0)
    positions = layout.pos[:, :2] * 0.9

    # picks = _check_channels_ordered(evoked.info, freqs)
    picks = _picks_to_idx(evoked.info, ['hbo', 'hbr'], exclude='bads', allow_empty=True)

    width, height = 0.05, 0.05

    for ii in picks:
        pos = positions[ii, :]

        # plot --- [lowerCorner_x, lowerCorner_y, width, height]
        ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])
        ax.plot(evoked.times, evoked.data[ii, :]*1e6, color=colour, linewidth=0.2)
        ax.fill_between(evoked.times, (evoked.data[ii, :]-error.data[ii, :])*1e6,
                        (evoked.data[ii, :]+error.data[ii, :])*1e6,
                        color=colour, alpha=0.1)
        ax.set_ylim(bottom=ylim[0], top=ylim[1])
        ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=3, length=2, width=0.5, labelcolor='w')
        ax.patch.set_alpha(0.3)
        ax.axhline(linewidth=0.1, color='k')
        # ax.set_title(f'{raw_long.ch_names[ii][:-4]}', fontsize=3)
        ax.set_title(f'{evoked.ch_names[ii][:-4]}', fontsize=3, pad=0)

    if add_legend is True:
        if hbdiff is True:
            custom_lines = [Line2D([0], [0], color=colour, lw=2)]
            fig.legend(custom_lines, ['ΔHbDiff'], loc='upper right')
        else:
            custom_lines = [Line2D([0], [0], color='r', lw=2),
                            Line2D([0], [0], color='b', lw=2)]
            fig.legend(custom_lines, ['ΔHbO', 'ΔHbR'], loc='upper right')

    # add an empty plot with labels
    ax = fig.add_axes([0.475, 0.075, 1 * width, 1 * height])
    ax.plot(evoked.times, np.zeros((len(evoked.times),)), '-w', linewidth=0.3)
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    ax.xaxis.set_major_locator(plt.MaxNLocator(2))
    ax.yaxis.set_major_locator(plt.MaxNLocator(2))
    ax.set_xlabel('Time (s)', fontsize=6)
    ax.set_ylabel('ΔHbO / ΔHbR (μM)', fontsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6)
    ax.axhline(linewidth=0.1, color='k')

    return fig


def plot_nirs_evoked(fig, evoked, colour, ylim, add_legend=False, hbdiff=False):
    """

    :param fig:
    :param evoked:
    :param colour:
    :param ylim:
    :param add_legend:
    :return:
    """

    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from mne.channels.layout import find_layout
    from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, \
        _check_channels_ordered
    import numpy as np
    from copy import deepcopy
    from mne.io.pick import _picks_to_idx

    layout = find_layout(evoked.info)
    layout = deepcopy(layout)
    layout.pos[:, :2] -= layout.pos[:, :2].min(0)
    layout.pos[:, :2] /= layout.pos[:, :2].max(0)
    positions = layout.pos[:, :2] * 0.9

    # picks = _check_channels_ordered(evoked.info, freqs)
    picks = _picks_to_idx(evoked.info, ['hbo', 'hbr'], exclude='bads', allow_empty=True)

    width, height = 0.05, 0.05

    for ii in picks:
        pos = positions[ii, :]

        # plot --- [lowerCorner_x, lowerCorner_y, width, height]
        ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])
        ax.plot(evoked.times, evoked.data[ii, :]*1e6, color=colour, linewidth=0.2)
        ax.set_ylim(bottom=ylim[0], top=ylim[1])
        ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=3, length=2, width=0.5, labelcolor='w')
        ax.patch.set_alpha(0.1)
        ax.axhline(linewidth=0.1, color='k')
        # ax.set_title(f'{raw_long.ch_names[ii][:-4]}', fontsize=3)
        ax.set_title(f'{evoked.ch_names[ii][:-4]}', fontsize=3, pad=0)

    if add_legend is True:
        if hbdiff is True:
            custom_lines = [Line2D([0], [0], color=colour, lw=2)]
            fig.legend(custom_lines, ['ΔHbDiff'], loc='lower right')
        else:
            custom_lines = [Line2D([0], [0], color='r', lw=2),
                            Line2D([0], [0], color='b', lw=2)]
            fig.legend(custom_lines, ['ΔHbO', 'ΔHbR'], loc='lower right')

    # add an empty plot with labels
    ax = fig.add_axes([0.075, 0.075, 1.5 * width, 1.5 * height])
    ax.plot(evoked.times, np.zeros((len(evoked.times),)), '-w', linewidth=0.3)
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    ax.xaxis.set_major_locator(plt.MaxNLocator(2))
    ax.yaxis.set_major_locator(plt.MaxNLocator(2))
    ax.set_xlabel('Time (s)', fontsize=6)
    ax.set_ylabel('ΔHbO / ΔHbR (μM)', fontsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6)
    ax.axhline(linewidth=0.1, color='k')

    return fig


def plot_nirs_evoked_error_mismatch(fig, evoked, error, colour, ylim, ch_idxs, add_legend=False, hbdiff=False):
    """

    :param fig:
    :param evoked:
    :param error:
    :param colour:
    :param ylim:
    :param add_legend:
    :return:
    """

    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from mne.channels.layout import find_layout
    from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, \
        _check_channels_ordered
    import numpy as np
    from copy import deepcopy
    from mne.io.pick import _picks_to_idx

    layout = find_layout(evoked.info)

    layout = deepcopy(layout)
    layout.pos[:, :2] -= layout.pos[:, :2].min(0)
    layout.pos[:, :2] /= layout.pos[:, :2].max(0)
    positions = layout.pos[:, :2] * 0.9

    # picks = _check_channels_ordered(evoked.info, freqs)
    picks_all = _picks_to_idx(evoked.info, ['hbo', 'hbr'], exclude='bads', allow_empty=True)

    width, height = 0.05, 0.05

    for ii, ch_idx in enumerate(ch_idxs):
        pos = positions[ch_idx, :]

        # plot --- [lowerCorner_x, lowerCorner_y, width, height]
        ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])
        ax.plot(evoked.times, evoked.data[ii, :]*1e6, color=colour, linewidth=0.2)
        ax.fill_between(evoked.times, (evoked.data[ii, :]-error.data[ii, :])*1e6,
                        (evoked.data[ii, :]+error.data[ii, :])*1e6,
                        color=colour, alpha=0.1)
        ax.set_ylim(bottom=ylim[0], top=ylim[1])
        ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=3, length=2, width=0.5, labelcolor='w')
        ax.patch.set_alpha(0.3)
        ax.axhline(linewidth=0.1, color='k')
        # ax.set_title(f'{raw_long.ch_names[ii][:-4]}', fontsize=3)
        ax.set_title(f'{evoked.ch_names[ch_idx][:-4]}', fontsize=3, pad=0)

    if add_legend is True:
        if hbdiff is True:
            custom_lines = [Line2D([0], [0], color=colour, lw=2)]
            fig.legend(custom_lines, ['ΔHbDiff'], loc='upper right')
        else:
            custom_lines = [Line2D([0], [0], color='r', lw=2),
                            Line2D([0], [0], color='b', lw=2)]
            fig.legend(custom_lines, ['ΔHbO', 'ΔHbR'], loc='upper right')

    # add an empty plot with labels
    ax = fig.add_axes([0.475, 0.075, 1 * width, 1 * height])
    ax.plot(evoked.times, np.zeros((len(evoked.times),)), '-w', linewidth=0.3)
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    ax.xaxis.set_major_locator(plt.MaxNLocator(2))
    ax.yaxis.set_major_locator(plt.MaxNLocator(2))
    ax.set_xlabel('Time (s)', fontsize=6)
    ax.set_ylabel('ΔHbO / ΔHbR (μM)', fontsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6)
    ax.axhline(linewidth=0.1, color='k')

    return fig
