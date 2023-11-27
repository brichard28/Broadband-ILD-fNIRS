# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 09:47:09 2023

@author: benri
"""

# Code to analyze SRM-NIRS-EEG-1
# Compares block averages and runs GLM on data from SRM-NIRS-EEG-1

# Import common libraries
import numpy as np
import pandas as pd
from itertools import compress
from collections import defaultdict
from copy import deepcopy
from pprint import pprint

# Import MNE processing
from mne.viz import plot_compare_evokeds
from mne.preprocessing.nirs import optical_density, beer_lambert_law, temporal_derivative_distribution_repair, scalp_coupling_index
from mne.io import read_raw_nirx

# Import MNE-NIRS processing
from mne_nirs.statistics import run_glm
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import statsmodels_to_results
from mne_nirs.channels import get_short_channels, get_long_channels
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.visualisation import plot_glm_group_topo
from mne_nirs.datasets import fnirs_motor_group
from mne_nirs.visualisation import plot_glm_surface_projection
from mne_nirs.io.fold import fold_channel_specificity
from mne_nirs.signal_enhancement import short_channel_regression

# Import MNE-BIDS processing
from mne_bids import BIDSPath, read_raw_bids, get_entity_vals

# Import StatsModels
import statsmodels.formula.api as smf

# Import Plotting Library
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scipy
import pandas

from itertools import compress

import mne
import os

import obspy
from obspy.signal.detrend import polynomial



root = 'C:/Users/benri/Nextcloud/nirs/data/'
all_fnirs_data_folders = ['/home/ben/Nextcloud/data/nirs/data/2023-09-21/2023-09-21_001',
'/home/ben/Nextcloud/data/nirs/data/2023-09-25/2023-09-25_001',
'/home/ben/Nextcloud/data/nirs/data/2023-09-26/2023-09-26_001',
'/home/ben/Nextcloud/data/nirs/data/2023-09-26/2023-09-26_002',
'/home/ben/Nextcloud/data/nirs/data/2023-10-02/2023-10-02_002',
'/home/ben/Nextcloud/data/nirs/data/2023-10-03/2023-10-03_001',
'/home/ben/Nextcloud/data/nirs/data/2023-10-03/2023-10-03_002',
'/home/ben/Nextcloud/data/nirs/data/2023-10-05/2023-10-05_001',
'/home/ben/Nextcloud/data/nirs/data/2023-10-10/2023-10-10_001',
'/home/ben/Nextcloud/data/nirs/data/2023-10-17/2023-10-17_001',
'/home/ben/Nextcloud/data/nirs/data/2023-10-19/2023-10-19_001',
'/home/ben/Nextcloud/data/nirs/data/2023-10-19/2023-10-19_002',
'/home/ben/Nextcloud/data/nirs/data/2023-10-26/2023-10-26_001',
'/home/ben/Nextcloud/data/nirs/data/2023-11-06/2023-11-06_001',
'/home/ben/Nextcloud/data/nirs/data/2023-11-09/2023-11-09_001',
'/home/ben/Nextcloud/data/nirs/data/2023-11-13/2023-11-13_001',
'/home/ben/Nextcloud/data/nirs/data/2023-11-16/2023-11-16_001']


# Before Control Condition: 'NDARYZ656HJ9','NDARCD778KPR','NDARMY829TKN','NDARLU426TBZ',
# Before Control Condition: '/home/ben/Nextcloud/data/nirs/data/2023-07-20/2023-07-20_002',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-07-20/2023-07-20_003',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-07-26/2023-07-26_001',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-08-02/2023-08-02_001',

# After ITD500 Control Condition added
subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARWK546QR2','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6',
                'NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3']
curr_subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARWK546QR2','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6',
                'NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3']


def individual_analysis(fnirs_data_folder, ID):
    
    os.environ["OMP_NUM_THREADS"] = "1"

    raw_intensity = read_raw_nirx(fnirs_data_folder, verbose=False, preload=True)
    
    # Before Control Condition
    #raw_intensity.annotations.rename({'15.0':'Masker Noise ITD 50 Targ Left',
    #                                  '14.0':'Masker Noise ILD 10 Targ Right',
    #                                  '13.0':'Masker Noise ILD 10 Targ Left',
    #                                  '12.0':'Masker Speech ILD 10 Targ Right',
    #                                  '11.0':'Masker Speech ITD 500 Targ Left',
    #                                  '10.0':'Masker Noise ITD 500 Targ Left',
    #                                  '9.0':'Masker Speech ITD 50 Targ Left',
    #                                  '8.0':'Masker Noise ITD 500 Targ Right',
    #                                  '7.0':'Masker Speech ITD 50 Targ Right',
    #                                  '6.0':'Masker Speech ILD 10 Targ Left',
    #                                  '5.0':'Masker Speech ITD 500 Targ Right',
    #                                  '4.0':'Masker Noise ITD 50 Targ Right',
    #                                  '3.0':'Hold',
    #                                  '2.0':'Exhale',
    #                                  '1.0':'Inhale'})
    raw_intensity.annotations.rename({'1.0':'Inhale',
                                      '2.0':'Exhale',
                                      '3.0':'Hold',
                                      '4.0':'m_speech__ild_0__itd_500__targ_r__control_0',
                                      '5.0':'m_noise__ild_0__itd_50__targ_l__control_0',
                                      '6.0':'m_noise__ild_0__itd_50__targ_r__control_0',
                                      '7.0':'m_speech__ild_70n__itd_0__targ_r__control_0',
                                      '8.0':'m_speech__ild_0__itd_50__targ_l__control_0',
                                      '9.0':'m_speech__ild_10__itd_0__targ_l__control_0',
                                      '10.0':'m_speech__ild_0__itd_500__targ_l__control_0',
                                      '11.0':'m_speech__ild_0__itd_500__targ_l__control_1',
                                      '12.0':'m_noise__ild_0__itd_500__targ_r__control_1',
                                      '13.0':'m_noise__ild_0__itd_500__targ_r__control_0',
                                      '14.0':'m_noise__ild_70n__itd_0__targ_l__control_0',
                                      '15.0':'m_noise__ild_10__itd_0__targ_l__control_0',
                                      '16.0':'m_speech__ild_10__itd_0__targ_r__control_0',
                                      '17.0':'m_noise__ild_10__itd_0__targ_r__control_0',
                                      '18.0':'m_speech__ild_0__itd_50__targ_r__control_0',
                                       '19.0':'m_speech__ild_0__itd_500__targ_r__control_1',
                                      '20.0':'m_noise__ild_70n__itd_0__targ_r__control_0',
                                      '21.0':'m_noise__ild_0__itd_500__targ_l__control_1',
                                      '22.0':'m_noise__ild_0__itd_500__targ_l__control_0',
                                      '23.0':'m_speech__ild_70n__itd_0__targ_l__control_0'})

    raw_od = optical_density(raw_intensity)
    
    # Scalp Coupling Index, label bad channels
    sci = scalp_coupling_index(raw_od)

    raw_od.info['bads'] = list(compress(raw_od.ch_names,sci < 0.8))

   # raw_od = short_channel_regression(raw_od, max_dist=0.01)
    raw_od = temporal_derivative_distribution_repair(raw_od, verbose=False)

    raw_od.interpolate_bads(verbose=False)

    raw_od.resample(10)

    raw_haemo = beer_lambert_law(raw_od, ppf=6)
    
    # Cut out just the short channels for creating a GLM repressor
    sht_chans = get_short_channels(raw_haemo)
    raw_haemo = get_long_channels(raw_haemo)
    
    raw_haemo._data = scipy.signal.detrend(raw_haemo._data, axis=-1, type='linear')

    # Filter data
    iir_params = dict({"order":3,"ftype":"butter","padlen":10000})
    raw_haemo = raw_haemo.filter(0.03, 0.7, iir_params=iir_params, method='iir', verbose=False)
    
    # # Plot events
    events, event_dict = mne.events_from_annotations(raw_haemo, verbose=False)
    # mne.viz.plot_events(events, event_id=event_dict,sfreq=raw_haemo.info['sfreq'])
    
       #reject_criteria = dict(hbo=10.0e-6)
       #flat_criteria = dict(hbo=2e-7)
    tmin, tmax = -5, 25
       
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject_by_annotation=True,
                           proj=True, baseline=(-5, 1), preload=False,
                           detrend=None, verbose=False)


    # APPLY BREATH HOLD NORMALIZATION
    
    # View Consistency across channels
    #fig2, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 6))
    #clims = dict(hbo=[-10, 10], hbr=[-10, 10])
    #epochs["m_noise__ild_0__itd_500__targ_l__control_1"].average().plot_image(axes=axes[:, 0], clim=clims)
    #epochs["m_noise__ild_0__itd_500__targ_r__control_1"].average().plot_image(axes=axes[:, 1], clim=clims)
    #for column, condition in enumerate(["m_noise__ild_0__itd_500__targ_l__control_1", "m_noise__ild_0__itd_500__targ_r__control_1"]):
    #    for ax in axes[:, column]:
    #        ax.set_title("{}: {}".format(condition, ax.get_title()))
            
    # Save Block Averages
    epochs_this_subject = epochs.to_data_frame()
    epochs_this_subject.to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " block averages.csv")

    # Create a design matrix
    design_matrix = make_first_level_design_matrix(raw_haemo, stim_dur=15.0, hrf_model='glover', drift_model='cosine'); #, drift_model= 'cosine')

    # Exclude bad epochs from design matrix
    
    # Append short channels mean to design matrix
    design_matrix["ShortHbO"] = np.mean(sht_chans.copy().pick(picks="hbo").get_data(), axis=0)
    design_matrix["ShortHbR"] = np.mean(sht_chans.copy().pick(picks="hbr").get_data(), axis=0)

    # Run GLM
    glm_est = run_glm(raw_haemo, design_matrix, noise_model='ar1',verbose=0)

    # Extract channel metrics
    cha = glm_est.to_dataframe()

    # Compute region of interest results from channel data
    #roi = glm_est.to_dataframe_region_of_interest(groups,
    #                                              design_matrix.columns,
    #                                              demographic_info=True)

    # Define left vs right contrast
    contrast_matrix = np.eye(design_matrix.shape[1])
    basic_conts = dict([(column, contrast_matrix[i])
                        for i, column in enumerate(design_matrix.columns)])
    contrast_LvR = basic_conts['m_speech__ild_0__itd_500__targ_l__control_1'] - basic_conts['m_noise__ild_0__itd_500__targ_l__control_1']

    # Compute defined contrast
    contrast = glm_est.compute_contrast(contrast_LvR)
    con = contrast.to_dataframe()

    # Add the participant ID to the dataframes
    cha["ID"] = con["ID"] = ID

    # Convert to uM for nicer plotting below.
    cha["theta"] = [t * 1.e6 for t in cha["theta"]]
    #roi["theta"] = [t * 1.e6 for t in roi["theta"]]
    con["effect"] = [t * 1.e6 for t in con["effect"]]
    
    return raw_haemo, cha, con, epochs



## Run analysis on all participants
df_roi = pd.DataFrame()  # To store region of interest results
df_cha = pd.DataFrame()  # To store channel level results
df_con = pd.DataFrame()  # To store channel level contrast results
all_hold_max = []

all_evokeds = defaultdict(list)    
    

for sub in curr_subject_ID:
    print("Subject:")  
    print(sub)

    # Create path to file based on experiment info
    which_folder = subject_ID.index(sub)
    fnirs_data_folder = all_fnirs_data_folders[which_folder]
    print(fnirs_data_folder)

    # Analyse data and return both ROI and channel results
    raw_haemo, channel, con, epochs = individual_analysis(fnirs_data_folder, sub)

    # Append individual results to all participants
    #df_roi = pd.concat([df_roi, roi], ignore_index=True)
    df_cha = pd.concat([df_cha, channel], ignore_index=True)
    df_con = pd.concat([df_con, con], ignore_index=True)




# Group Results
grp_results = df_cha.query("Condition in ['Inhale','Exhale','Hold','m_speech__ild_0__itd_500__targ_r__control_0','m_noise__ild_0__itd_50__targ_l__control_0','m_noise__ild_0__itd_50__targ_r__control_0','m_speech__ild_70n__itd_0__targ_r__control_0','m_speech__ild_0__itd_50__targ_l__control_0','m_speech__ild_10__itd_0__targ_l__control_0','m_speech__ild_0__itd_500__targ_l__control_0','m_speech__ild_0__itd_500__targ_l__control_1','m_noise__ild_0__itd_500__targ_r__control_1','m_noise__ild_0__itd_500__targ_r__control_0','m_noise__ild_70n__itd_0__targ_l__control_0','m_noise__ild_10__itd_0__targ_l__control_0','m_speech__ild_10__itd_0__targ_r__control_0','m_noise__ild_10__itd_0__targ_r__control_0','m_speech__ild_0__itd_50__targ_r__control_0','m_speech__ild_0__itd_500__targ_r__control_1','m_noise__ild_70n__itd_0__targ_r__control_0','m_noise__ild_0__itd_500__targ_l__control_1','m_noise__ild_0__itd_500__targ_l__control_0','m_speech__ild_70n__itd_0__targ_l__control_0']")
grp_results = grp_results.query("Chroma in ['hbo']")


roi_model = smf.mixedlm("theta ~ Condition:ID",
                        grp_results, groups=grp_results["ID"]).fit(method='nm')
print(roi_model.summary())

grp_results.to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1.csv")


