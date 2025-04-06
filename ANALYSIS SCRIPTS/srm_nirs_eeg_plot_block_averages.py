# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_block_averages(hbo_speech_masker_filename ,hbr_speech_masker_filename, hbo_noise_masker_filename ,hbr_noise_masker_filename):
    # Use a breakpoint in the code line below to debug your script.
    hbo_speech_masker_csv = pd.read_csv(hbo_speech_masker_filename)
    hbr_speech_masker_csv = pd.read_csv(hbr_speech_masker_filename)
    hbo_noise_masker_csv = pd.read_csv(hbo_noise_masker_filename)
    hbr_noise_masker_csv = pd.read_csv(hbr_noise_masker_filename)

    all_hbo_data_pfc_speech_masker_small_itd = []
    all_hbo_data_pfc_speech_masker_large_itd = []
    all_hbo_data_pfc_speech_masker_natural_ild = []
    all_hbo_data_pfc_speech_masker_broadband_ild = []

    all_hbo_data_stg_speech_masker_small_itd = []
    all_hbo_data_stg_speech_masker_large_itd = []
    all_hbo_data_stg_speech_masker_natural_ild = []
    all_hbo_data_stg_speech_masker_broadband_ild = []

    all_hbr_data_pfc_speech_masker_small_itd = []
    all_hbr_data_pfc_speech_masker_large_itd = []
    all_hbr_data_pfc_speech_masker_natural_ild = []
    all_hbr_data_pfc_speech_masker_broadband_ild = []

    all_hbr_data_stg_speech_masker_small_itd = []
    all_hbr_data_stg_speech_masker_large_itd = []
    all_hbr_data_stg_speech_masker_natural_ild = []
    all_hbr_data_stg_speech_masker_broadband_ild = []

    all_hbo_data_pfc_noise_masker_small_itd = []
    all_hbo_data_pfc_noise_masker_large_itd = []
    all_hbo_data_pfc_noise_masker_natural_ild = []
    all_hbo_data_pfc_noise_masker_broadband_ild = []

    all_hbo_data_stg_noise_masker_small_itd = []
    all_hbo_data_stg_noise_masker_large_itd = []
    all_hbo_data_stg_noise_masker_natural_ild = []
    all_hbo_data_stg_noise_masker_broadband_ild = []

    all_hbr_data_pfc_noise_masker_small_itd = []
    all_hbr_data_pfc_noise_masker_large_itd = []
    all_hbr_data_pfc_noise_masker_natural_ild = []
    all_hbr_data_pfc_noise_masker_broadband_ild = []

    all_hbr_data_stg_noise_masker_small_itd = []
    all_hbr_data_stg_noise_masker_large_itd = []
    all_hbr_data_stg_noise_masker_natural_ild = []
    all_hbr_data_stg_noise_masker_broadband_ild = []

    curr_subjects = np.unique(hbo_speech_masker_csv['S'])
    curr_channels = np.unique(hbo_speech_masker_csv['Channel'])
    for isub, sub in enumerate(curr_subjects):
        for ichannel, channel in enumerate(curr_channels):
            this_hbo_data_speech_masker = hbo_speech_masker_csv.loc[
                (hbo_speech_masker_csv['S'] == sub) & (hbo_speech_masker_csv['Channel'] == channel)]
            this_hbr_data_speech_masker = hbr_speech_masker_csv.loc[
                (hbr_speech_masker_csv['S'] == sub) & (hbr_speech_masker_csv['Channel'] == channel)]
            this_hbo_data_noise_masker = hbo_noise_masker_csv.loc[
                (hbo_noise_masker_csv['S'] == sub) & (hbo_noise_masker_csv['Channel'] == channel)]
            this_hbr_data_noise_masker = hbr_noise_masker_csv.loc[
                (hbr_noise_masker_csv['S'] == sub) & (hbr_noise_masker_csv['Channel'] == channel)]
            if ichannel > 6:
                all_hbo_data_stg_speech_masker_small_itd.append(this_hbo_data_speech_masker['0'].to_numpy())
                all_hbo_data_stg_speech_masker_large_itd.append(this_hbo_data_speech_masker['1'].to_numpy())
                all_hbo_data_stg_speech_masker_natural_ild.append(this_hbo_data_speech_masker['2'].to_numpy())
                all_hbo_data_stg_speech_masker_broadband_ild.append(this_hbo_data_speech_masker['3'].to_numpy())

                all_hbr_data_stg_speech_masker_small_itd.append(this_hbr_data_speech_masker['0'].to_numpy())
                all_hbr_data_stg_speech_masker_large_itd.append(this_hbr_data_speech_masker['1'].to_numpy())
                all_hbr_data_stg_speech_masker_natural_ild.append(this_hbr_data_speech_masker['2'].to_numpy())
                all_hbr_data_stg_speech_masker_broadband_ild.append(this_hbr_data_speech_masker['3'].to_numpy())

                all_hbo_data_stg_noise_masker_small_itd.append(this_hbo_data_noise_masker['0'].to_numpy())
                all_hbo_data_stg_noise_masker_large_itd.append(this_hbo_data_noise_masker['1'].to_numpy())
                all_hbo_data_stg_noise_masker_natural_ild.append(this_hbo_data_noise_masker['2'].to_numpy())
                all_hbo_data_stg_noise_masker_broadband_ild.append(this_hbo_data_noise_masker['3'].to_numpy())

                all_hbr_data_stg_noise_masker_small_itd.append(this_hbr_data_noise_masker['0'].to_numpy())
                all_hbr_data_stg_noise_masker_large_itd.append(this_hbr_data_noise_masker['1'].to_numpy())
                all_hbr_data_stg_noise_masker_natural_ild.append(this_hbr_data_noise_masker['2'].to_numpy())
                all_hbr_data_stg_noise_masker_broadband_ild.append(this_hbr_data_noise_masker['3'].to_numpy())

            elif ichannel <= 6:
                all_hbo_data_pfc_speech_masker_small_itd.append(this_hbo_data_speech_masker['0'].to_numpy())
                all_hbo_data_pfc_speech_masker_large_itd.append(this_hbo_data_speech_masker['1'].to_numpy())
                all_hbo_data_pfc_speech_masker_natural_ild.append(this_hbo_data_speech_masker['2'].to_numpy())
                all_hbo_data_pfc_speech_masker_broadband_ild.append(this_hbo_data_speech_masker['3'].to_numpy())

                all_hbr_data_pfc_speech_masker_small_itd.append(this_hbr_data_speech_masker['0'].to_numpy())
                all_hbr_data_pfc_speech_masker_large_itd.append(this_hbr_data_speech_masker['1'].to_numpy())
                all_hbr_data_pfc_speech_masker_natural_ild.append(this_hbr_data_speech_masker['2'].to_numpy())
                all_hbr_data_pfc_speech_masker_broadband_ild.append(this_hbr_data_speech_masker['3'].to_numpy())

                all_hbo_data_pfc_noise_masker_small_itd.append(this_hbo_data_noise_masker['0'].to_numpy())
                all_hbo_data_pfc_noise_masker_large_itd.append(this_hbo_data_noise_masker['1'].to_numpy())
                all_hbo_data_pfc_noise_masker_natural_ild.append(this_hbo_data_noise_masker['2'].to_numpy())
                all_hbo_data_pfc_noise_masker_broadband_ild.append(this_hbo_data_noise_masker['3'].to_numpy())

                all_hbr_data_pfc_noise_masker_small_itd.append(this_hbr_data_noise_masker['0'].to_numpy())
                all_hbr_data_pfc_noise_masker_large_itd.append(this_hbr_data_noise_masker['1'].to_numpy())
                all_hbr_data_pfc_noise_masker_natural_ild.append(this_hbr_data_noise_masker['2'].to_numpy())
                all_hbr_data_pfc_noise_masker_broadband_ild.append(this_hbr_data_noise_masker['3'].to_numpy())


    epoch_start_time = -5
    epoch_end_time = 20
    num_timepoints = 255 #(np.abs(epoch_start_time) + np.abs(epoch_end_time))*10.2
    single_block_time = np.linspace(epoch_start_time,epoch_end_time,num_timepoints)


    fig, block_average_axes = plt.subplots(nrows = 2, ncols = 4)

    plot_indices_order = [[0,0],[0,1],[1,0],[1,1],[0,2],[0,3],[1,2],[1,3]]
    iplot = 0
    for masker_type in ['speech','noise']:
        for roi in ['pfc','stg']:
            for chroma in ['hbo','hbr']:
                this_color = 'r' if chroma == 'hbo' else 'b'
                this_row = 0 if masker_type == 'speech' else 1
                if roi == 'pfc':
                    these_column_choices = [0,1]
                else:
                    these_column_choices = [2,3]

                this_data = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_small_itd"), axis=0)
                this_error = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_small_itd"), axis=0)/(np.sqrt(len(curr_subjects)) - 1)
                block_average_axes[this_row,these_column_choices[0]].plot(single_block_time, this_data, linestyle='-', color=this_color)
                block_average_axes[this_row,these_column_choices[0]].fill_between(single_block_time, this_data + this_error, this_data - this_error,
                                                      alpha=0.1, color=this_color)

                this_data = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_large_itd"), axis=0)
                this_error = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_large_itd"), axis=0) / (
                            np.sqrt(len(curr_subjects)) - 1)
                block_average_axes[this_row, these_column_choices[0]].plot(single_block_time, this_data, linestyle='--',
                                                                           color=this_color)
                block_average_axes[this_row, these_column_choices[0]].fill_between(single_block_time,
                                                                                   this_data + this_error,
                                                                                   this_data - this_error,
                                                                                   alpha=0.1, color=this_color)

                this_data = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_natural_ild"), axis=0)
                this_error = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_natural_ild"), axis=0) / (
                            np.sqrt(len(curr_subjects)) - 1)
                block_average_axes[this_row, these_column_choices[1]].plot(single_block_time, this_data, linestyle='-',
                                                                           color=this_color)
                block_average_axes[this_row, these_column_choices[1]].fill_between(single_block_time,
                                                                                   this_data + this_error,
                                                                                   this_data - this_error,
                                                                                   alpha=0.1, color=this_color)

                this_data = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_broadband_ild"), axis=0)
                this_error = np.nanmean(eval(f"all_{chroma}_data_{roi}_{masker_type}_masker_broadband_ild"), axis=0) / (
                        np.sqrt(len(curr_subjects)) - 1)
                block_average_axes[this_row, these_column_choices[1]].plot(single_block_time, this_data, linestyle='--',
                                                                           color=this_color)
                block_average_axes[this_row, these_column_choices[1]].fill_between(single_block_time,
                                                                                   this_data + this_error,
                                                                                   this_data - this_error,
                                                                                   alpha=0.1, color=this_color)



    for ax in block_average_axes.flat:  # Iterates through all subplots
        ax.set_ylim([-0.05, 0.15])
    plt.savefig("/home/ben/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/python_block_averages.png")



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    hbo_speech_masker_filename = "/home/ben/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/all_subjects_uncorr_block_average_speech_masker.csv"
    hbr_speech_masker_filename = "/home/ben/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/all_subjects_uncorr_block_average_hbr_speech_masker.csv"
    hbo_noise_masker_filename = "/home/ben/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/all_subjects_uncorr_block_average_noise_masker.csv"
    hbr_noise_masker_filename = "/home/ben/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/all_subjects_uncorr_block_average_hbr_noise_masker.csv"
    plot_block_averages(hbo_speech_masker_filename ,hbr_speech_masker_filename, hbo_noise_masker_filename ,hbr_noise_masker_filename)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
