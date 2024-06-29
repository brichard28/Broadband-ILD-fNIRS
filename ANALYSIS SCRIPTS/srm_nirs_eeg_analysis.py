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

import numpy as np
import mne
import pickle
import matplotlib
# matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import os
import pandas as pd
from collections import defaultdict

mne.set_config('MNE_BROWSER_BACKEND', 'qt')
#from nirx_movement import mark_aux_movement_bad
import mne_nirs
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
from mne_nirs.channels import (get_long_channels,
                               get_short_channels)
#import pandas as pd
#from nilearn.plotting import plot_design_matrix
from run_preproc_NIRS import preprocess_NIRX
from mne_nirs.io.snirf import read_snirf_aux_data
#from sklearn.preprocessing import StandardScaler, MinMaxScaler
#from scipy import signal
from scipy import stats

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
data_root + '2024-04-16/2024-04-16_001',
data_root + '2024-04-24/2024-04-24_001',
data_root + '2024-04-24/2024-04-24_002',
data_root + '2024-04-25/2024-04-25_001',
data_root + '2024-04-25/2024-04-25_002']


# Before Control Condition: 'NDARYZ656HJ9','NDARCD778KPR','NDARMY829TKN','NDARLU426TBZ',
# Before Control Condition: '/home/ben/Nextcloud/data/nirs/data/2023-07-20/2023-07-20_002',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-07-20/2023-07-20_003',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-07-26/2023-07-26_001',
                            #'/home/ben/Nextcloud/data/nirs/data/2023-08-02/2023-08-02_001',

# After ITD500 Control Condition added
# All subject IDs
subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARFV472HU7',
                'NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARSZ622LR8','NDARGS283RM9','NDARRED356WS','NDARHUG535MO','NDARFIN244AL','NDARKAI888JU','NDARGRA444CE','NDARBAA679HA']
# The subjects we would like to run right now
curr_subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARFV472HU7','NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARGS283RM9','NDARRED356WS', 'NDARHUG535MO','NDARFIN244AL','NDARKAI888JU','NDARBAA679HA']

masker_type = 'speech' # type of masker to analyze on this run
glm_dur = 7

n_subjects = len(all_fnirs_data_folders)

n_long_channels = 14
fs = 10.2
n_timepoints = 306
task_type = 'Ben_SvN'

def individual_analysis(fnirs_data_folder, ID):
    
    os.environ["OMP_NUM_THREADS"] = "1"

    base_dir = "C:/Users/elibu/Documents/NIRx/Data/"
    save_dir = "C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/Data/"
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    plot_dir = "C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/Plots/"
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)
    
    # Read in the raw nirs file
    #raw_intensity = read_raw_nirx(fnirs_data_folder, verbose=False, preload=True)
    data = mne.io.read_raw_nirx(f"{fnirs_data_folder}/{fnirs_data_folder[-14:]}_config.hdr",
                                verbose=False, preload=True)

    data_snirf = mne.io.read_raw_snirf(f"{fnirs_data_folder}/{fnirs_data_folder[-14:]}.snirf",
                                       optode_frame="mri", preload=True)

    aux_snirf = read_snirf_aux_data(f"{fnirs_data_folder}/{fnirs_data_folder[-14:]}.snirf",
                                    data_snirf)
    
    # ---------------------------------------------------------------
    # -----------------      Preprocess the Data            ---------
    # ---------------------------------------------------------------  
    if masker_type == 'speech': # speech masker
        data.annotations.rename({'1.0':'Inhale',
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
        data_snirf.annotations.rename({'1':'Inhale',
                                          '2':'Exhale',
                                          '3':'Hold',
                                          '4':'speech',
                                          '5':'ild_0__itd_50',
                                          '6':'ild_0__itd_50',
                                          '7':'speech',
                                          '8':'speech',
                                          '9':'speech',
                                          '10':'speech',
                                          '11':'control',
                                          '12':'control',
                                          '13':'ild_0__itd_500',
                                          '14':'ild_70n__itd_0',
                                          '15':'ild_10__itd_0',
                                          '16':'speech',
                                          '17':'ild_10__itd_0',
                                          '18':'speech',
                                          '19':'control',
                                          '20':'ild_70n__itd_0',
                                          '21':'control',
                                          '22':'ild_0__itd_500',
                                          '23':'speech'})
    elif masker_type == 'noise': # noise masker
        data.annotations.rename({'1.0':'Inhale',
                                          '2.0':'Exhale',
                                          '3.0':'Hold',
                                          '4.0':'ild_0__itd_500',
                                          '5.0':'noise',
                                          '6.0':'noise',
                                          '7.0':'ild_70n__itd_0',
                                          '8.0':'ild_0__itd_50',
                                          '9.0':'ild_10__itd_0',
                                          '10.0':'ild_0__itd_500',
                                          '11.0':'control',
                                          '12.0':'control',
                                          '13.0':'noise',
                                          '14.0':'noise',
                                          '15.0':'noise',
                                          '16.0':'ild_10__itd_0',
                                          '17.0':'noise',
                                          '18.0':'ild_0__itd_50',
                                          '19.0':'control',
                                          '20.0':'noise',
                                          '21.0':'control',
                                          '22.0':'noise',
                                          '23.0':'ild_70n__itd_0'})
        data_snirf.annotations.rename({'1':'Inhale',
                                          '2':'Exhale',
                                          '3':'Hold',
                                          '4':'ild_0__itd_500',
                                          '5':'noise',
                                          '6':'noise',
                                          '7':'ild_70n__itd_0',
                                          '8':'ild_0__itd_50',
                                          '9':'ild_10__itd_0',
                                          '10':'ild_0__itd_500',
                                          '11':'control',
                                          '12':'control',
                                          '13':'noise',
                                          '14':'noise',
                                          '15':'noise',
                                          '16':'ild_10__itd_0',
                                          '17':'noise',
                                          '18':'ild_0__itd_50',
                                          '19':'control',
                                          '20':'noise',
                                          '21':'control',
                                          '22':'noise',
                                          '23':'ild_70n__itd_0'})
    else: # both maskers:
        data.annotations.rename({'1.0': 'Inhale',
                                  '2.0': 'Exhale',
                                  '3.0': 'Hold',
                                  '4.0': 'ild_0__itd_500',
                                  '5.0': 'ild_0__itd_50',
                                '6.0': 'ild_0__itd_50',
                                  '7.0': 'ild_70n__itd_0',
                                  '8.0': 'ild_0__itd_50',
                                  '9.0': 'ild_10__itd_0',
                                  '10.0': 'ild_0__itd_500',
                                  '11.0': 'control',
                                  '12.0': 'control',
                                  '13.0': 'ild_0__itd_500',
                                  '14.0': 'ild_70n__itd_0',
                                  '15.0': 'ild_10__itd_0',
                                  '16.0': 'ild_10__itd_0',
                                  '17.0': 'ild_10__itd_0',
                                  '18.0': 'ild_0__itd_50',
                                  '19.0': 'control',
                                  '20.0': 'ild_70n__itd_0',
                                  '21.0': 'control',
                                  '22.0': 'ild_0__itd_500',
                                  '23.0': 'ild_70n__itd_0'})

        data_snirf.annotations.rename({'1': 'Inhale',
                                        '2': 'Exhale',
                                        '3': 'Hold',
                                        '4': 'ild_0__itd_500',
                                        '5': 'ild_0__itd_50',
                                        '6': 'ild_0__itd_50',
                                        '7': 'ild_70n__itd_0',
                                        '8': 'ild_0__itd_50',
                                        '9': 'ild_10__itd_0',
                                        '10': 'ild_0__itd_500',
                                        '11': 'control',
                                        '12': 'control',
                                        '13': 'ild_0__itd_500',
                                        '14': 'ild_70n__itd_0',
                                        '15': 'ild_10__itd_0',
                                        '16': 'ild_10__itd_0',
                                        '17': 'ild_10__itd_0',
                                        '18': 'ild_0__itd_50',
                                        '19': 'control',
                                        '20': 'ild_70n__itd_0',
                                        '21': 'control',
                                        '22': 'ild_0__itd_500',
                                        '23': 'ild_70n__itd_0'})

   #  # Convert to optical density
   #  raw_od = optical_density(raw_intensity)
    
   #  # Scalp Coupling Index, label bad channels
   #  sci = scalp_coupling_index(raw_od)

   #  # Add 'bads' to info
    
   #  raw_od.info['bads'] = list(compress(raw_od.ch_names,sci < 0.8))

   # # raw_od = short_channel_regression(raw_od, max_dist=0.01)
   
   #  # Apply TDDR
   #  raw_od = temporal_derivative_distribution_repair(raw_od, verbose=False)
    
   #  # Resample to 3 Hz
   #  raw_od.resample(3) # 10
    
   #  # Create separate object for block averages (will run short channel on these, but use short channels as a regressor in the GLM for betas)
   #  raw_od_for_block_averages = short_channel_regression(raw_od, max_dist=0.01)
    
   #  #raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
   #  #raw_haemo_for_block_averages = beer_lambert_law(raw_od_for_block_averages, ppf=0.1)
   #  raw_haemo = mne_modified_beer_lambert_law(raw_od) # TRYING ELIS FUNCTION
   #  raw_haemo_for_block_averages = mne_modified_beer_lambert_law(raw_od_for_block_averages)
    
   #  # Cut out just the short channels for creating a GLM repressor
   #  sht_chans = get_short_channels(raw_haemo)
   #  raw_haemo = get_long_channels(raw_haemo)
   #  raw_haemo_for_block_averages = get_long_channels(raw_haemo_for_block_averages)
    
   #  bad_channels = list(compress(raw_haemo.ch_names,sci < 0.8))
    
   #  # Filter data
   #  iir_params = dict({"order":3,"ftype":"butter","padlen":10000})
   #  raw_haemo = raw_haemo.filter(0.01, 0.3, iir_params=iir_params, method='iir', verbose=False)
   #  raw_haemo_for_block_averages = raw_haemo_for_block_averages.filter(0.01, 0.3, iir_params=iir_params, method='iir', verbose=False)
   #  #raw_haemo.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005) # 0.01, 0.3, 0.2, 0.005
   #  #raw_haemo_for_block_averages.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005)
       
   #  if ID == 'NDARBA306US5' or ID == 'NDARDC882NK4':
   #      events_to_remove = [45,46,47]
   #      raw_haemo.annotations.delete(events_to_remove)
   #      raw_haemo_for_block_averages.annotations.delete(events_to_remove)
   #  if ID == 'NDARAZC45TW3':
   #      events_to_remove = [1,2]
   #      raw_haemo.annotations.delete(events_to_remove)
   #      raw_haemo_for_block_averages.annotations.delete(events_to_remove)
        
    tmin, tmax = -5, 25 # timings for block averages
    
    
    events, event_dict = mne.events_from_annotations(data, verbose=False)
    
    

    raw_haemo, null = preprocess_NIRX(data, data_snirf, event_dict,
                                           save=True,
                                           savename=save_dir + f'{ID}_{task_type}_preproc_nirs.fif',
                                           plot_steps=False,
                                           crop=False, crop_low=0, crop_high=0,
                                           events_modification=False, reject=True,
                                           short_regression=True, events_from_snirf=False,
                                           drop_short=False, negative_enhancement=False,
                                           snr_thres=3, filter_type='iir')
    raw_haemo_for_block_averages = raw_haemo.copy()
    
    sht_chans = get_short_channels(raw_haemo)
    raw_haemo = get_long_channels(raw_haemo)
    raw_haemo_for_block_averages = get_long_channels(raw_haemo_for_block_averages)
    
    # Define epochs
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject_by_annotation=True,
                           proj=True, baseline=(-5, 0), preload=True,
                           detrend=None, verbose=False)
    # APPLY BREATH HOLD NORMALIZATION
    # Get data from first inhale to end of hold
    hold_length = 30
    epochs_for_hold = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=0, tmax=hold_length,reject_by_annotation=True,
                           proj=True, baseline=None, preload=True,
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
    #print(f"Max During Breathing = {max_during_breathing}")
    
    # Redefine epochs
    reject_criteria = dict(hbo=5e-6) # 10
    #flat_criteria = dict(hbo=0.1e-6)
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject=reject_criteria,reject_by_annotation=True,
                           proj=True, baseline=(-5, 0), preload=True,
                           detrend=None, verbose=False)

    # Remove rejected epochs from design matrix
    # epochs_to_remove = []
    # peak_to_peak_threshold = 10
    # for iepoch in range(len(epochs.drop_log)):
    #     z = np.max(epochs[iepoch].get_data(),axis=2) - np.min(epochs[iepoch].get_data(),axis=2) > peak_to_peak_threshold #z = np.std(epochs[iepoch].get_data(),axis=2) > 3  #z = np.max(epochs[iepoch].get_data(),axis=2) > 2
    #     if np.size(epochs.drop_log[iepoch]) > 0:
    #         epochs_to_remove.append(iepoch)
    #         print(f"deleted epoch: {iepoch}")
        
    #     elif z.any():
    #         epochs_to_remove.append(iepoch)
    #         print(f"deleted epoch: {iepoch}")
    # raw_haemo.annotations.delete(epochs_to_remove)
    # raw_haemo_for_block_averages.annotations.delete(epochs_to_remove)
    
    # Save Block Averages
    epochs_for_block_averages = mne.Epochs(raw_haemo_for_block_averages, events, event_id=event_dict,
                           tmin=tmin, tmax=tmax,reject=reject_criteria,reject_by_annotation=True,
                           proj=True, baseline=(-5, 0), preload=True,
                           detrend=None, verbose=False)

    #epochs_for_block_averages.drop_channels(bad_channels)
    
    epochs_this_subject = epochs_for_block_averages.to_data_frame()
    
    # PLOTTING 
    # epochs_for_block_averages["speech"].plot_image(combine="mean",vmin=-1,vmax=1,ts_args=dict(ylim=dict(hbo=[-1, 1], hbr=[-1, 1])))
    # epochs_for_block_averages["noise"].plot_image(combine="mean",vmin=-1,vmax=1,ts_args=dict(ylim=dict(hbo=[-1, 1], hbr=[-1, 1])))

    
    if user == 'Laptop':
        epochs_this_subject.to_csv("C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " block averages SPEECH NO BREATH.csv")
        #pandas.DataFrame(epochs_to_remove).to_csv("C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " epochs NOISE NO BREATH.csv")
    else:
        epochs_this_subject.to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " block averages SPEECH NO BREATH.csv")
       # pandas.DataFrame(epochs_to_remove).to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/" + sub + " epochs NOISE NO BREATH.csv")

    
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
    
    grp_results.to_csv("C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 collapsed attend and masker PFC time constant SPEECH NO BREATH.csv")
else:

    grp_results.to_csv("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 collapsed attend and masker PFC time constant SPEECH NO BREATH.csv")


