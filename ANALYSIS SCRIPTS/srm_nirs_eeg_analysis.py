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

from mne_modified_beer_lambert_law import mne_modified_beer_lambert_law

from sklearn.preprocessing import StandardScaler, MinMaxScaler

# Define Subject Files
root = ''
user = 'Laptop'
if user == 'Laptop':
    data_root = 'C:/Users/benri/Downloads/'

else:
    data_root = '/home/ben/Nextcloud/data/nirs/data/'
    
all_fnirs_data_folders = [data_root + '2023-09-21/2023-09-21_001',
data_root + '2023-09-25/2023-09-25_001',
data_root + '2023-09-26/2023-09-26_001',
data_root + '2023-09-26/2023-09-26_002',
data_root + '2023-10-02/2023-10-02_002',
data_root + '2023-10-03/2023-10-03_001',
data_root + '2023-10-05/2023-10-05_001',
data_root + '2023-10-10/2023-10-10_001',
data_root + '2023-10-17/2023-10-17_001',
data_root + '2023-10-19/2023-10-19_001',
data_root + '2023-10-19/2023-10-19_002',
data_root + '2023-10-24/2023-10-24_001',
data_root + '2023-10-26/2023-10-26_001',
data_root + '2023-11-06/2023-11-06_001',
data_root + '2023-11-09/2023-11-09_001',
data_root + '2023-11-13/2023-11-13_001',
data_root + '2023-11-16/2023-11-16_001',
data_root + '2023-11-20/2023-11-20_001',
data_root + '2023-11-27/2023-11-27_001',
data_root + '2023-12-04/2023-12-04_001',
data_root + '2024-01-30/2024-01-30_001',
data_root + '2024-01-31/2024-01-31_001',
data_root + '2024-02-01/2024-02-01_001',
data_root + '2024-04-09/2024-04-09_001',
data_root + '2024-04-16/2024-04-16_001']


# Before Control Condition: 'NDARYZ656HJ9','NDARCD778KPR','NDARMY829TKN','NDARLU426TBZ',
# Before Control Condition: '/home/ben/Nextcloud/data/nirs/data/2023-07-20/2023-07-20_002',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-07-20/2023-07-20_003',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-07-26/2023-07-26_001',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-08-02/2023-08-02_001',

# After ITD500 Control Condition added
# All subject IDs
subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARFV472HU7',
                'NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARSZ622LR8','NDARGS283RM9','NDARRED356WS','NDARHUG535MO']
# The subjects we would like to run right now
curr_subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARFV472HU7','NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARGS283RM9','NDARRED356WS', 'NDARHUG535MO']


def individual_analysis(fnirs_data_folder, ID):
    
    os.environ["OMP_NUM_THREADS"] = "1"

    # Read in the raw nirs file
    raw_intensity = read_raw_nirx(fnirs_data_folder, verbose=False, preload=True)
    
    # Before Control Condition
    #raw_intensity.annotations.rename({'1.0':'Inhale',
                                      # '2.0':'Exhale',
                                      # '3.0':'Hold',
                                      # '4.0':'m_speech__ild_0__itd_500__targ_r__control_0',
                                      # '5.0':'m_noise__ild_0__itd_50__targ_l__control_0',
                                      # '6.0':'m_noise__ild_0__itd_50__targ_r__control_0',
                                      # '7.0':'m_speech__ild_70n__itd_0__targ_r__control_0',
                                      # '8.0':'m_speech__ild_0__itd_50__targ_l__control_0',
                                      # '9.0':'m_speech__ild_10__itd_0__targ_l__control_0',
                                      # '10.0':'m_speech__ild_0__itd_500__targ_l__control_0',
                                      # '11.0':'m_speech__ild_0__itd_500__targ_l__control_1',
                                      # '12.0':'m_noise__ild_0__itd_500__targ_r__control_1',
                                      # '13.0':'m_noise__ild_0__itd_500__targ_r__control_0',
                                      # '14.0':'m_noise__ild_70n__itd_0__targ_l__control_0',
                                      # '15.0':'m_noise__ild_10__itd_0__targ_l__control_0',
                                      # '16.0':'m_speech__ild_10__itd_0__targ_r__control_0',
                                      # '17.0':'m_noise__ild_10__itd_0__targ_r__control_0',
                                      # '18.0':'m_speech__ild_0__itd_50__targ_r__control_0',
                                      #  '19.0':'m_speech__ild_0__itd_500__targ_r__control_1',
                                      # '20.0':'m_noise__ild_70n__itd_0__targ_r__control_0',
                                      # '21.0':'m_noise__ild_0__itd_500__targ_l__control_1',
                                      # '22.0':'m_noise__ild_0__itd_500__targ_l__control_0',
                                      # '23.0':'m_speech__ild_70n__itd_0__targ_l__control_0'})
    # raw_intensity.annotations.rename({'1.0':'Inhale',
    #                                   '2.0':'Exhale',
    #                                   '3.0':'Hold',
    #                                   '4.0':'m_speech__ild_0__itd_500__targ_r__control_0',
    #                                   '5.0':'m_noise__ild_0__itd_50__targ_l__control_0',
    #                                   '6.0':'m_noise__ild_0__itd_50__targ_r__control_0',
    #                                   '7.0':'m_speech__ild_70n__itd_0__targ_r__control_0',
    #                                   '8.0':'m_speech__ild_0__itd_50__targ_l__control_0',
    #                                   '9.0':'m_speech__ild_10__itd_0__targ_l__control_0',
    #                                   '10.0':'m_speech__ild_0__itd_500__targ_l__control_0',
    #                                   '11.0':'m_speech__ild_0__itd_500__targ_l__control_1',
    #                                   '12.0':'m_noise__ild_0__itd_500__targ_r__control_1',
    #                                   '13.0':'m_noise__ild_0__itd_500__targ_r__control_0',
    #                                   '14.0':'m_noise__ild_70n__itd_0__targ_l__control_0',
    #                                   '15.0':'m_noise__ild_10__itd_0__targ_l__control_0',
    #                                   '16.0':'m_speech__ild_10__itd_0__targ_r__control_0',
    #                                   '17.0':'m_noise__ild_10__itd_0__targ_r__control_0',
    #                                   '18.0':'m_speech__ild_0__itd_50__targ_r__control_0',
    #                                    '19.0':'m_speech__ild_0__itd_500__targ_r__control_1',
    #                                   '20.0':'m_noise__ild_70n__itd_0__targ_r__control_0',
    #                                   '21.0':'m_noise__ild_0__itd_500__targ_l__control_1',
    #                                   '22.0':'m_noise__ild_0__itd_500__targ_l__control_0',
    #                                   '23.0':'m_speech__ild_70n__itd_0__targ_l__control_0'})
    
    # CHANGED HERE TO COLLAPSE ACROSS ATTEND LEFT AND RIGHT, SPEECH VS. NOISE FOR PFC PURPOSES
    # raw_intensity.annotations.rename({'1.0':'Inhale',
    #                                   '2.0':'Exhale',
    #                                   '3.0':'Hold',
    #                                   '4.0':'ild_0__itd_500',
    #                                   '5.0':'ild_0__itd_50',
    #                                   '6.0':'ild_0__itd_50',
    #                                   '7.0':'ild_70n__itd_0',
    #                                   '8.0':'ild_0__itd_50',
    #                                   '9.0':'ild_10__itd_0',
    #                                   '10.0':'ild_0__itd_500',
    #                                   '11.0':'control',
    #                                   '12.0':'control',
    #                                   '13.0':'ild_0__itd_500',
    #                                   '14.0':'ild_70n__itd_0',
    #                                   '15.0':'ild_10__itd_0',
    #                                   '16.0':'ild_10__itd_0',
    #                                   '17.0':'ild_10__itd_0',
    #                                   '18.0':'ild_0__itd_50',
    #                                     '19.0':'control',
    #                                   '20.0':'ild_70n__itd_0',
    #                                   '21.0':'control',
    #                                   '22.0':'ild_0__itd_500',
    #                                   '23.0':'ild_70n__itd_0'})
    
    # raw_intensity.annotations.rename({'1.0':'Inhale',
    #                                   '2.0':'Exhale',
    #                                   '3.0':'Hold',
    #                                   '4.0':'ild_0__itd_500',
    #                                   '5.0':'noise',
    #                                   '6.0':'noise',
    #                                   '7.0':'ild_70n__itd_0',
    #                                   '8.0':'ild_0__itd_50',
    #                                   '9.0':'ild_10__itd_0',
    #                                   '10.0':'ild_0__itd_500',
    #                                   '11.0':'control',
    #                                   '12.0':'control',
    #                                   '13.0':'noise',
    #                                   '14.0':'noise',
    #                                   '15.0':'noise',
    #                                   '16.0':'ild_10__itd_0',
    #                                   '17.0':'noise',
    #                                   '18.0':'ild_0__itd_50',
    #                                   '19.0':'control',
    #                                   '20.0':'noise',
    #                                   '21.0':'control',
    #                                   '22.0':'noise',
    #                                   '23.0':'ild_70n__itd_0'})
    
    raw_intensity.annotations.rename({'1.0':'Inhale',
                                      '2.0':'Exhale',
                                      '3.0':'Hold',
                                      '4.0':'speech',
                                      '5.0':'ild_0__itd_50',
                                      '6.0':'ild_0__itd_50',
                                      '7.0':'speech',
                                      '8.0':'speech',
                                      '9.0':'speech',
                                      '10.0':'speech',
                                      '11.0':'control',
                                      '12.0':'control',
                                      '13.0':'ild_0__itd_500',
                                      '14.0':'ild_70n__itd_0',
                                      '15.0':'ild_10__itd_0',
                                      '16.0':'speech',
                                      '17.0':'ild_10__itd_0',
                                      '18.0':'speech',
                                        '19.0':'control',
                                      '20.0':'ild_70n__itd_0',
                                      '21.0':'control',
                                      '22.0':'ild_0__itd_500',
                                      '23.0':'speech'})
    
    # raw_intensity.annotations.rename({'1.0':'Inhale',
    #                                     '2.0':'Exhale',
    #                                     '3.0':'Hold',
    #                                     '4.0':'speech',
    #                                     '5.0':'noise',
    #                                     '6.0':'noise',
    #                                     '7.0':'speech',
    #                                     '8.0':'speech',
    #                                     '9.0':'speech',
    #                                     '10.0':'speech',
    #                                     '11.0':'control',
    #                                     '12.0':'control',
    #                                     '13.0':'noise',
    #                                     '14.0':'noise',
    #                                     '15.0':'noise',
    #                                     '16.0':'speech',
    #                                     '17.0':'noise',
    #                                     '18.0':'speech',
    #                                     '19.0':'control',
    #                                     '20.0':'noise',
    #                                     '21.0':'control',
    #                                     '22.0':'noise',
    #                                     '23.0':'speech'})

    # Convert to optical density
    raw_od = optical_density(raw_intensity)
    
    # Scalp Coupling Index, label bad channels
    sci = scalp_coupling_index(raw_od)

    # Add 'bads' to info
    raw_od.info['bads'] = list(compress(raw_od.ch_names,sci < 0.8))

   # raw_od = short_channel_regression(raw_od, max_dist=0.01)
   
    # Apply TDDR
    raw_od = temporal_derivative_distribution_repair(raw_od, verbose=False)
    
    # # Interpolate bad channels
    raw_od.interpolate_bads(verbose=False)

    # Resample to 3 Hz
    raw_od.resample(3) # 10
    
    # Create separate object for block averages (will run short channel on these, but use short channels as a regressor in the GLM for betas)
    raw_od_for_block_averages = short_channel_regression(raw_od, max_dist=0.01)
    
    #raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
    #raw_haemo_for_block_averages = beer_lambert_law(raw_od_for_block_averages, ppf=0.1)
    raw_haemo = mne_modified_beer_lambert_law(raw_od) # TRYING ELIS FUNCTION
    raw_haemo_for_block_averages = mne_modified_beer_lambert_law(raw_od_for_block_averages)
    
    # Cut out just the short channels for creating a GLM repressor
    sht_chans = get_short_channels(raw_haemo)
    raw_haemo = get_long_channels(raw_haemo)
    raw_haemo_for_block_averages = get_long_channels(raw_haemo_for_block_averages)
    
    # Filter data
    #iir_params = dict({"order":3,"ftype":"butter","padlen":10000})
    #raw_haemo = raw_haemo.filter(0.03, 0.7, iir_params=iir_params, method='iir', verbose=False)
    raw_haemo.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005) # 0.01, 0.3, 0.2, 0.005
    raw_haemo_for_block_averages.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005)
       
    if ID == 'NDARBA306US5' or ID == 'NDARDC882NK4':
        events_to_remove = [45,46,47]
        raw_haemo.annotations.delete(events_to_remove)
        raw_haemo_for_block_averages.annotations.delete(events_to_remove)
    if ID == 'NDARAZC45TW3':
        events_to_remove = [1,2]
        raw_haemo.annotations.delete(events_to_remove)
        raw_haemo_for_block_averages.annotations.delete(events_to_remove)
        
    tmin, tmax = -5, 35 # timings for block averages
    
    # Shift events by 1.5 seconds (because of cue)
    #raw_haemo.annotations.onset += 1.5
    
    
    # Redefine events
    events, event_dict = mne.events_from_annotations(raw_haemo, verbose=False)
    
    # Define epochs
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject_by_annotation=True,
                           proj=True, baseline=(-5, 0), preload=False,
                           detrend=None, verbose=False)
    # APPLY BREATH HOLD NORMALIZATION
    # Get data from first inhale to end of hold
    hold_length = 30
    epochs_for_hold = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=0, tmax=hold_length,reject_by_annotation=True,
                           proj=True, baseline=None, preload=False,
                           detrend=None, verbose=False)
    hold_data = epochs_for_hold["Hold"].get_data(picks='hbo')
    hold_data = np.mean(hold_data,axis=(0))
    max_during_breathing = np.max(hold_data[:,45:91], axis=1) * 1e6
    max_during_breathing = np.repeat(max_during_breathing,2,axis=0)
    z = raw_haemo.get_data() 
    zz = raw_haemo_for_block_averages.get_data()
    #raw_haemo._data = np.divide(z,max_during_breathing[:, None])
    #raw_haemo_for_block_averages._data = np.divide(zz,max_during_breathing[:, None])
    # BREATH HOLD IS NOT CHANNEL SPECIFIC RIGHT NOW NEEDS TO CHANGE
    print(f"Max During Breathing = {max_during_breathing}")
    
    # Redefine epochs
    reject_criteria = dict(hbo=2e6)
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject=reject_criteria,reject_by_annotation=True,
                           proj=True, baseline=(-5, 0), preload=False,
                           detrend=None, verbose=False)

    # Remove rejected epochs from design matrix
    epochs_to_remove = []
    peak_to_peak_threshold = 10
    for iepoch in range(len(epochs.drop_log)):
        z = np.max(epochs[iepoch].get_data(),axis=2) - np.min(epochs[iepoch].get_data(),axis=2) > peak_to_peak_threshold #z = np.std(epochs[iepoch].get_data(),axis=2) > 3  #z = np.max(epochs[iepoch].get_data(),axis=2) > 2
        if np.size(epochs.drop_log[iepoch]) > 0:
            epochs_to_remove.append(iepoch)
            print(f"deleted epoch: {iepoch}")
        
        elif z.any():
            epochs_to_remove.append(iepoch)
            print(f"deleted epoch: {iepoch}")
    raw_haemo.annotations.delete(epochs_to_remove)
    raw_haemo_for_block_averages.annotations.delete(epochs_to_remove)
    
    # Save Block Averages
    epochs_for_block_averages = mne.Epochs(raw_haemo_for_block_averages, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject=reject_criteria,reject_by_annotation=True,
                           proj=True, baseline=(-5, 0), preload=False,
                           detrend=None, verbose=False)
    epochs_this_subject = epochs_for_block_averages.to_data_frame()
    
    if user == 'Laptop':
        epochs_this_subject.to_csv("C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " block averages NOISE NO BREATH.csv")
        pandas.DataFrame(epochs_to_remove).to_csv("C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " epochs NOISE NO BREATH.csv")
    else:
        epochs_this_subject.to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " block averages NOISE NO BREATH.csv")
        pandas.DataFrame(epochs_to_remove).to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " epochs NOISE NO BREATH.csv")

    
    # Create a design matrix
    design_matrix_hbo = make_first_level_design_matrix(raw_haemo.pick(picks='hbo'),drift_model=None,high_pass=0.01,hrf_model='spm',stim_dur=10.8)

    design_matrix_hbo["Linear"] = np.arange(0, np.shape(design_matrix_hbo)[0])
    
    # Append short channels mean to design matrix
    design_matrix_hbo["ShortHbO"] = np.mean(sht_chans.copy().pick(picks="hbo").get_data(),axis=0)
    
    #design_matrix["ShortHbR"] = np.mean(sht_chans.copy().pick(picks="hbr").get_data(), axis=0)

    #min_max_scaler = MinMaxScaler()
    #X_minmax = min_max_scaler.fit_transform(design_matrix_hbo)
    #design_matrix_min_max = pd.DataFrame(X_minmax, columns=design_matrix_hbo.columns.tolist())
    
    
    # Run GLM
    raw_haemo_glm = raw_haemo.copy().pick(picks='hbo')
    glm_est = run_glm(raw_haemo_glm, design_matrix_hbo, noise_model='ar1',verbose=0)

    # Extract channel metrics
    cha = glm_est.to_dataframe()

    # Add the participant ID to the dataframes
    cha["ID"] = ID

    # Convert to uM for nicer plotting below.
    cha["theta"] = [t * 1.e6 for t in cha["theta"]]
    #roi["theta"] = [t * 1.e6 for t in roi["theta"]]
    #con["effect"] = [t * 1.e6 for t in con["effect"]]
    
    return raw_haemo, cha, epochs



## Run analysis on all participants
df_roi = pd.DataFrame()  # To store region of interest results
df_cha = pd.DataFrame()  # To store channel level results
all_hold_max = []

all_evokeds = defaultdict(list)    
    

for sub in curr_subject_ID: # For each subject...
    # Print subject ID
    print("Subject:")  
    print(sub)

    # Create path to file based on experiment info
    which_folder = subject_ID.index(sub)
    fnirs_data_folder = all_fnirs_data_folders[which_folder]
    print(fnirs_data_folder)

    # Analyse data and return both ROI and channel results
    raw_haemo, channel, epochs = individual_analysis(fnirs_data_folder, sub)

    # Append individual results to all participants
    df_cha = pd.concat([df_cha, channel], ignore_index=True)




# Group Results
#grp_results = df_cha.query("Condition in ['Inhale','Exhale','Hold','m_speech__ild_0__itd_500__targ_r__control_0','m_noise__ild_0__itd_50__targ_l__control_0','m_noise__ild_0__itd_50__targ_r__control_0','m_speech__ild_70n__itd_0__targ_r__control_0','m_speech__ild_0__itd_50__targ_l__control_0','m_speech__ild_10__itd_0__targ_l__control_0','m_speech__ild_0__itd_500__targ_l__control_0','m_speech__ild_0__itd_500__targ_l__control_1','m_noise__ild_0__itd_500__targ_r__control_1','m_noise__ild_0__itd_500__targ_r__control_0','m_noise__ild_70n__itd_0__targ_l__control_0','m_noise__ild_10__itd_0__targ_l__control_0','m_speech__ild_10__itd_0__targ_r__control_0','m_noise__ild_10__itd_0__targ_r__control_0','m_speech__ild_0__itd_50__targ_r__control_0','m_speech__ild_0__itd_500__targ_r__control_1','m_noise__ild_70n__itd_0__targ_r__control_0','m_noise__ild_0__itd_500__targ_l__control_1','m_noise__ild_0__itd_500__targ_l__control_0','m_speech__ild_70n__itd_0__targ_l__control_0']")
grp_results = df_cha.query("Condition in ['Inhale','Exhale','Hold','control','ild_0__itd_50','ild_0__itd_500','ild_10__itd_0','ild_70n__itd_0']")
#grp_results = df_cha.query("Condition in ['Inhale','Exhale','Hold','control','speech','noise']")
grp_results = grp_results.query("Chroma in ['hbo']")

if user == 'Laptop':
    
    grp_results.to_csv("C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 collapsed attend and masker PFC time constant NOISE NO BREATH.csv")
else:

    grp_results.to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 collapsed attend and masker PFC time constant NOISE NO BREATH.csv")


