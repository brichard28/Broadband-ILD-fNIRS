# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:21:52 2024

@author: Benjamin Richardson
"""

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


# Ben preprocessing specific imports
from mne.preprocessing.nirs import optical_density, beer_lambert_law, temporal_derivative_distribution_repair, scalp_coupling_index
from itertools import compress
from mne_nirs.signal_enhancement import short_channel_regression
from mne_modified_beer_lambert_law import mne_modified_beer_lambert_law
from mne.io import read_raw_nirx
from mpl_toolkits.mplot3d import Axes3D
# ---------------------------------------------------------------
# -----------------          Data Parameters            ---------
# ---------------------------------------------------------------
wdir = os.path.dirname(__file__)

# Define Subject Files
data_root = 'C:/Users/benri/Downloads/' #'C:/Users/elibu/Documents/NIRx/Data/Ben_SvN/'

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
                          data_root + '2023-10-26/2023-10-26_001',
                          data_root + '2023-11-06/2023-11-06_001',
                          data_root + '2023-11-09/2023-11-09_001',
                          data_root + '2023-11-13/2023-11-13_001',
                          data_root + '2023-11-16/2023-11-16_001',
                          data_root + '2023-11-20/2023-11-20_001',
                          data_root + '2023-11-27/2023-11-27_001',
                          data_root + '2023-12-04/2023-12-04_001',
                          data_root + '2024-01-30/2024-01-30_001',
                          data_root + '2024-02-01/2024-02-01_001',
                          data_root + '2024-04-09/2024-04-09_001',
                          data_root + '2024-04-16/2024-04-16_001',
                          data_root + '2024-04-24/2024-04-24_001',
                          data_root + '2024-04-24/2024-04-24_002',
                          data_root + '2024-04-25/2024-04-25_002',
                          data_root + '2024-06-14/2024-06-14_001',
                          data_root + '2024-07-01/2024-07-01_002',
                          data_root + '2024-07-02/2024-07-02_001',
                          data_root + '2024-07-02/2024-07-02_002']

subject_ID = ['NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5',
                'NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARDC882NK4',
                'NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3',
                'NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8',
                'NDARLJ581GD7','NDARGS283RM9','NDARRED356WS',
                'NDARHUG535MO',
                'NDARFIN244AL','NDARKAI888JU','NDARBAA679HA',
                'NDARUXL573SS',
                'NDARMOL966PB','NDARGHM426BL','NDARSEW256ZA']  

target_direction = 'left' # type of masker to analyze on this run
masker_type = 'noise'
glm_dur = 7
preprocessing_type = "Ben"

n_subjects = len(all_fnirs_data_folders)
n_long_channels = 14
fs = 10.2
n_timepoints = 255



# set up the arrays to hold all subject data
subject_data_hold = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd50 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild70n = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild10 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild70n_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild10_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild70n_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild10_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild70n_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild10_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)


subject_data_hold_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd50_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild70n_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild10_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_itd500_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild70n_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild10_GLM = np.full((n_subjects, n_long_channels), np.nan)



# put into a larger array with all subjects data!
subject_data_itd50_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_itd500_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild70n_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild10_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)

num_channels_removed = np.full(n_subjects, np.nan)
age = np.full(n_subjects, np.nan)
sex = np.full(n_subjects, np.nan)

range_BH_response = np.zeros((n_subjects, n_long_channels))

all_evokeds = defaultdict(list)
subject_info = []
# loop through all subjects and all sessions (takes a while)
for ii, subject_num in enumerate(range(n_subjects)):
    
    
    
    os.environ["OMP_NUM_THREADS"] = "1"
    
    subject = subject_ID[ii]
    task_type = 'Ben_SvN'
    base_dir = "C:/Users/elibu/Documents/NIRx/Data/"
    save_dir = "C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/Data/"
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    plot_dir = "C:/Users/benri/Documents/GitHub/SRM-NIRS-EEG/ANALYSIS SCRIPTS/Eli Analysis/Plots/"
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)

    plot_steps = False

    # ---------------------------------------------------------------
    # -----------------      Load the Data        ---------
    # ---------------------------------------------------------------
    if preprocessing_type == "Eli":
        data = mne.io.read_raw_nirx(f"{all_fnirs_data_folders[ii]}/{all_fnirs_data_folders[ii][-14:]}_config.hdr",
                                    verbose=False, preload=True)
    
        data_snirf = mne.io.read_raw_snirf(f"{all_fnirs_data_folders[ii]}/{all_fnirs_data_folders[ii][-14:]}.snirf",
                                           optode_frame="mri", preload=True)
    
        aux_snirf = read_snirf_aux_data(f"{all_fnirs_data_folders[ii]}/{all_fnirs_data_folders[ii][-14:]}.snirf",
                                        data_snirf)
    elif preprocessing_type == "Ben":
        data = read_raw_nirx(f"{all_fnirs_data_folders[ii]}/", verbose=False, preload=True)
        data_snirf = mne.io.read_raw_snirf(f"{all_fnirs_data_folders[ii]}/{all_fnirs_data_folders[ii][-14:]}.snirf",
                                           optode_frame="mri", preload=True)
    
    
    # ---------------------------------------------------------------
    # -----------------      Preprocess the Data            ---------
    # ---------------------------------------------------------------  
    if target_direction == 'right' and masker_type == 'speech':
        data.annotations.rename({'1.0':'Inhale',
                                '2.0':'Exhale',
                                '3.0':'Hold',
                                '4.0':'noise',
                                '5.0':'attend_left',
                                '6.0':'ild_0__itd_50',
                                '7.0':'noise',
                                '8.0':'noise',
                                '9.0':'noise',
                                '10.0':'noise',
                                '11.0':'control',
                                '12.0':'control',
                                '13.0':'ild_0__itd_500',
                                '14.0':'attend_left',
                                '15.0':'attend_left',
                                '16.0':'noise',
                                '17.0':'ild_10__itd_0',
                                '18.0':'noise',
                                '19.0':'control',
                                '20.0':'ild_70n__itd_0',
                                '21.0':'control',
                                '22.0':'attend_left',
                                '23.0':'noise'})
        data_snirf.annotations.rename({'1':'Inhale',
                                    '2':'Exhale',
                                    '3':'Hold',
                                    '4':'noise',
                                    '5':'attend_left',
                                    '6':'ild_0__itd_50',
                                    '7':'noise',
                                    '8':'noise',
                                    '9':'noise',
                                    '10':'noise',
                                    '11':'control',
                                    '12':'control',
                                    '13':'ild_0__itd_500',
                                    '14':'attend_left',
                                    '15':'attend_left',
                                    '16':'noise',
                                    '17':'ild_10__itd_0',
                                    '18':'noise',
                                    '19':'control',
                                    '20':'ild_70n__itd_0',
                                    '21':'control',
                                    '22':'attend_left',
                                    '23':'noise'})
    elif target_direction == 'right' and masker_type == 'noise':
        data.annotations.rename({'1.0':'Inhale',
                                          '2.0':'Exhale',
                                          '3.0':'Hold',
                                          '4.0':'ild_0__itd_500',
                                          '5.0':'noise',
                                          '6.0':'noise',
                                          '7.0':'ild_70n__itd_0',
                                          '8.0':'attend_left',
                                          '9.0':'attend_left',
                                          '10.0':'attend_left',
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
                                          '23.0':'attend_left'})
        data_snirf.annotations.rename({'1':'Inhale',
                                          '2':'Exhale',
                                          '3':'Hold',
                                          '4':'ild_0__itd_500',
                                          '5':'noise',
                                          '6':'noise',
                                          '7':'ild_70n__itd_0',
                                          '8':'attend_left',
                                          '9':'attend_left',
                                          '10':'attend_left',
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
                                          '23':'attend_left'})

    elif target_direction == 'left' and masker_type == 'speech':
        data.annotations.rename({'1.0':'Inhale',
                                          '2.0':'Exhale',
                                          '3.0':'Hold',
                                          '4.0':'speech',
                                          '5.0':'ild_0__itd_50',
                                          '6.0':'attend_right',
                                          '7.0':'speech',
                                          '8.0':'speech',
                                          '9.0':'speech',
                                          '10.0':'speech',
                                          '11.0':'control',
                                          '12.0':'control',
                                          '13.0':'attend_right',
                                          '14.0':'ild_70n__itd_0',
                                          '15.0':'ild_10__itd_0',
                                          '16.0':'speech',
                                          '17.0':'attend_right',
                                          '18.0':'speech',
                                            '19.0':'control',
                                          '20.0':'attend_right',
                                          '21.0':'control',
                                          '22.0':'ild_0__itd_500',
                                          '23.0':'speech'})
        data_snirf.annotations.rename({'1':'Inhale',
                                          '2':'Exhale',
                                          '3':'Hold',
                                          '4':'speech',
                                          '5':'ild_0__itd_50',
                                          '6':'attend_right',
                                          '7':'speech',
                                          '8':'speech',
                                          '9':'speech',
                                          '10':'speech',
                                          '11':'control',
                                          '12':'control',
                                          '13':'attend_right',
                                          '14':'ild_70n__itd_0',
                                          '15':'ild_10__itd_0',
                                          '16':'speech',
                                          '17':'attend_right',
                                          '18':'speech',
                                          '19':'control',
                                          '20':'attend_right',
                                          '21':'control',
                                          '22':'ild_0__itd_500',
                                          '23':'speech'})

    elif target_direction == 'left' and masker_type == 'noise':
        data.annotations.rename({'1.0':'Inhale',
                                          '2.0':'Exhale',
                                          '3.0':'Hold',
                                          '4.0':'attend_right',
                                          '5.0':'noise',
                                          '6.0':'noise',
                                          '7.0':'attend_right',
                                          '8.0':'ild_0__itd_50',
                                          '9.0':'ild_10__itd_0',
                                          '10.0':'ild_0__itd_500',
                                          '11.0':'control',
                                          '12.0':'control',
                                          '13.0':'noise',
                                          '14.0':'noise',
                                          '15.0':'noise',
                                          '16.0':'attend_right',
                                          '17.0':'noise',
                                          '18.0':'attend_right',
                                          '19.0':'control',
                                          '20.0':'noise',
                                          '21.0':'control',
                                          '22.0':'noise',
                                          '23.0':'ild_70n__itd_0'})
        data_snirf.annotations.rename({'1':'Inhale',
                                          '2':'Exhale',
                                          '3':'Hold',
                                          '4':'attend_right',
                                          '5':'noise',
                                          '6':'noise',
                                          '7':'attend_right',
                                          '8':'ild_0__itd_50',
                                          '9':'ild_10__itd_0',
                                          '10':'ild_0__itd_500',
                                          '11':'control',
                                          '12':'control',
                                          '13':'noise',
                                          '14':'noise',
                                          '15':'noise',
                                          '16':'attend_right',
                                          '17':'noise',
                                          '18':'attend_right',
                                          '19':'control',
                                          '20':'noise',
                                          '21':'control',
                                          '22':'noise',
                                          '23':'ild_70n__itd_0'})
        
    # ---------------------------------------------------------------
    # -------------               Preprocessing             ---------
    # ---------------------------------------------------------------
    
    if preprocessing_type == "Eli":
        events, event_dict = mne.events_from_annotations(data, verbose=False)
    
        raw_haemo_temp, null = preprocess_NIRX(data, data_snirf, event_dict,
                                               save=True,
                                               savename=save_dir + f'{subject}_{task_type}_preproc_nirs.fif',
                                               plot_steps=False,
                                               crop=False, crop_low=0, crop_high=0,
                                               events_modification=False, reject=True,
                                               short_regression=True, events_from_snirf=False,
                                               drop_short=False, negative_enhancement=False,
                                               snr_thres=3, filter_type='iir')
    
        raw_haemo_short = get_short_channels(raw_haemo_temp)
        raw_haemo_filt = get_long_channels(raw_haemo_temp)
        
        # extra_regressors = aux_snirf.reset_index(drop=True)
        #
        # order = 4  # You can adjust the order as needed
        # [b, a] = signal.iirfilter(N=order, Wn=0.01 / (0.5 * raw_haemo_filt.info['sfreq']), btype='high', ftype='butter')
        # filtered_signals = extra_regressors.iloc[:, 1:].apply(lambda col: signal.filtfilt(b, a, col), axis=0)

    elif preprocessing_type == "Ben":
        events, event_dict = mne.events_from_annotations(data, verbose=False)
        del event_dict['Hold']
        del event_dict['Inhale']
        del event_dict['Exhale']
       # del event_dict['control']
        # Convert to optical density
        raw_od = optical_density(data_snirf)
         
        # Scalp Coupling Index, label bad channels
        sci = scalp_coupling_index(raw_od)

        # Add 'bads' to info
         
        raw_od.info['bads'] = list(compress(raw_od.ch_names,sci < 0.35))
        
        

        # raw_od = short_channel_regression(raw_od, max_dist=0.01)
        
        # Apply TDDR
        raw_od = temporal_derivative_distribution_repair(raw_od, verbose=False)
         
        # Resample to 3 Hz
        #raw_od.resample(3) # 10
         
        # Create separate object for block averages (will run short channel on these, but use short channels as a regressor in the GLM for betas)
        raw_od_regressed = short_channel_regression(raw_od.copy(), max_dist=0.01)
         
        #raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
        raw_haemo = mne_modified_beer_lambert_law(raw_od_regressed) # TRYING ELIS FUNCTION
         
        # Filter data
        iir_params = dict({"order":3,"ftype":"butter","padlen":10000}) # 3
        raw_haemo = raw_haemo.filter(0.01, 0.3, iir_params=iir_params, method='iir', verbose=False) #0.01,0.3
        #raw_haemo.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005) # 0.01, 0.3, 0.2, 0.005
        
        raw_haemo_short = get_short_channels(raw_haemo)
        raw_haemo_filt = get_long_channels(raw_haemo)
        
        
        num_channels_removed[ii] = len(list(raw_haemo_filt.info['bads']))/2
        age[ii] = 2024 - data.info['subject_info']['birthday'][0]
        sex[ii] = data.info['subject_info']['sex']


    # ---------------------------------------------------------------
    # -------------               Epoching                  ---------
    # ---------------------------------------------------------------
    reject_criteria = dict(hbo=5e-6)#5e-6
    #flat_criteria = dict(hbo=0.05e-6)
    tmin, tmax = -5, 20

    epochs = mne.Epochs(raw_haemo_filt, events,  # events_block,
                        event_id=event_dict,  # event_dict_total,
                        tmin=tmin, tmax=tmax,
                        baseline= None, # (-5, 0)
                        reject = reject_criteria,
                       # flat = flat_criteria,
                        preload=True, detrend=None, verbose=True,
                        on_missing='warn')
    #epochs.plot_drop_log()
    #plt.show()
    epochs.drop_bad()
    
    
    n_conditions = 5
    conditions = ['Hold','ild_0__itd_50','ild_0__itd_500','ild_70n__itd_0','ild_10__itd_0']


    # mark where the bad channels are
    chan_hbo = epochs.copy().pick('hbo').info['ch_names']
    chan_hbo_bad = list(epochs.copy().pick('hbo').info['bads'])

    chan_indices_bad = [i for i in range(len(chan_hbo)) if chan_hbo[i] in chan_hbo_bad]
    chan_indices_good = [i for i in range(len(chan_hbo)) if chan_hbo[i] not in chan_hbo_bad]

    # split the data into speech, noise, control,
    data_itd50 = epochs["ild_0__itd_50"].get_data(picks='hbo')
    data_itd500 = epochs["ild_0__itd_500"].get_data(picks='hbo')
    data_ild70n = epochs["ild_70n__itd_0"].get_data(picks='hbo')
    data_ild10 = epochs["ild_10__itd_0"].get_data(picks='hbo')
    
    data_itd50_hbr = epochs["ild_0__itd_50"].get_data(picks='hbr')
    data_itd500_hbr = epochs["ild_0__itd_500"].get_data(picks='hbr')
    data_ild70n_hbr = epochs["ild_70n__itd_0"].get_data(picks='hbr')
    data_ild10_hbr = epochs["ild_10__itd_0"].get_data(picks='hbr')
    
    # ---------------------------------------------------------------
    # -----------------    Baselining and Averaging         ---------
    # ---------------------------------------------------------------
    data_itd50_avg = np.full((len(chan_indices_good), n_timepoints), np.nan)
    data_itd500_avg = np.full((len(chan_indices_good), n_timepoints), np.nan)
    data_ild70n_avg = np.full((len(chan_indices_good), n_timepoints), np.nan)
    data_ild10_avg = np.full((len(chan_indices_good), n_timepoints), np.nan)

    data_itd50_avg_hbr= np.full((len(chan_indices_good), n_timepoints), np.nan)
    data_itd500_avg_hbr= np.full((len(chan_indices_good), n_timepoints), np.nan)
    data_ild70n_avg_hbr= np.full((len(chan_indices_good), n_timepoints), np.nan)
    data_ild10_avg_hbr= np.full((len(chan_indices_good), n_timepoints), np.nan)
    
    all_data_hbo_this_subject = np.concatenate((data_itd50, data_itd500, data_ild70n, data_ild10) , axis = 0)
    all_data_hbr_this_subject = np.concatenate((data_itd50_hbr, data_itd500_hbr, data_ild70n_hbr, data_ild10_hbr) , axis = 0)
    
    for ichannel in range(len(chan_indices_good)):
        
        data_itd50_avg[ichannel,:] = np.nanmean(data_itd50[:,ichannel,:]  - np.nanmean(data_itd50[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        data_itd500_avg[ichannel,:] = np.nanmean(data_itd500[:,ichannel,:]  - np.nanmean(data_itd500[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        data_ild70n_avg[ichannel,:] = np.nanmean(data_ild70n[:,ichannel,:]  -  np.nanmean(data_ild70n[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        data_ild10_avg[ichannel,:] = np.nanmean(data_ild10[:,ichannel,:]  -  np.nanmean(data_ild10[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        
        data_itd50_avg_hbr[ichannel,:] = np.nanmean(data_itd50_hbr[:,ichannel,:]  -  np.nanmean(data_itd50_hbr[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        data_itd500_avg_hbr[ichannel,:] = np.nanmean(data_itd500_hbr[:,ichannel,:]  -  np.nanmean(data_itd500_hbr[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        data_ild70n_avg_hbr[ichannel,:] = np.nanmean(data_ild70n_hbr[:,ichannel,:]  -  np.nanmean(data_ild70n_hbr[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
        data_ild10_avg_hbr[ichannel,:] = np.nanmean(data_ild10_hbr[:,ichannel,:]  -  np.nanmean(data_ild10_hbr[:,ichannel,0:int(5*fs)], axis = (0,1)), axis=0)
    
    # for ichannel in range(len(chan_indices_good)):
        
    #     data_itd50_avg[ichannel,:] = np.nanmean(data_itd50[:,ichannel,:], axis=0)
    #     data_itd500_avg[ichannel,:] = np.nanmean(data_itd500[:,ichannel,:], axis=0)
    #     data_ild70n_avg[ichannel,:] = np.nanmean(data_ild70n[:,ichannel,:], axis=0)
    #     data_ild10_avg[ichannel,:] = np.nanmean(data_ild10[:,ichannel,:], axis=0)
        
    #     data_itd50_avg_hbr[ichannel,:] = np.nanmean(data_itd50_hbr[:,ichannel,:], axis=0)
    #     data_itd500_avg_hbr[ichannel,:] = np.nanmean(data_itd500_hbr[:,ichannel,:], axis=0)
    #     data_ild70n_avg_hbr[ichannel,:] = np.nanmean(data_ild70n_hbr[:,ichannel,:], axis=0)
    #     data_ild10_avg_hbr[ichannel,:] = np.nanmean(data_ild10_hbr[:,ichannel,:], axis=0)
    
    # need to mark the indices where the good channels are!

    # put into a larger array with all subjects data!
    subject_data_itd50[ii, chan_indices_good, :] = 1e6*data_itd50_avg.copy()
    subject_data_itd500[ii, chan_indices_good, :] = 1e6*data_itd500_avg.copy()
    subject_data_ild70n[ii, chan_indices_good, :] = 1e6*data_ild70n_avg.copy()
    subject_data_ild10[ii, chan_indices_good, :] = 1e6*data_ild10_avg.copy()
  
    subject_data_itd50_hbr[ii, chan_indices_good, :] = 1e6*data_itd50_avg_hbr.copy()
    subject_data_itd500_hbr[ii, chan_indices_good, :] = 1e6*data_itd500_avg_hbr.copy()
    subject_data_ild70n_hbr[ii, chan_indices_good, :] = 1e6*data_ild70n_avg_hbr.copy()
    subject_data_ild10_hbr[ii, chan_indices_good, :] = 1e6*data_ild10_avg_hbr.copy()
    
    subject_data_itd50_baselined = subject_data_itd50.copy()
    subject_data_itd500_baselined = subject_data_itd500.copy()
    subject_data_ild70n_baselined = subject_data_ild70n.copy()
    subject_data_ild10_baselined = subject_data_ild10.copy()
    
    subject_data_itd50_hbr_baselined = subject_data_itd50_hbr.copy()
    subject_data_itd500_hbr_baselined = subject_data_itd500_hbr.copy()
    subject_data_ild70n_hbr_baselined = subject_data_ild70n_hbr.copy()
    subject_data_ild10_hbr_baselined = subject_data_ild10_hbr.copy()
    

    # for ichannel in range(n_long_channels):
    #     subject_data_itd50_baselined[ii,ichannel,:] = (subject_data_itd50[ii,ichannel,:]  - np.nanmean(subject_data_itd50[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_itd50[:,:,0:int(5*fs)], axis=(0,1,2))
    #     subject_data_itd500_baselined[ii,ichannel,:]  = (subject_data_itd500[ii,ichannel,:]  - np.nanmean(subject_data_itd500[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_itd500[:,:,0:int(5*fs)], axis=(0,1,2))
    #     subject_data_ild70n_baselined[ii,ichannel,:]  = (subject_data_ild70n[ii,ichannel,:]  - np.nanmean(subject_data_ild70n[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_ild70n[:,:,0:int(5*fs)], axis=(0,1,2))
    #     subject_data_ild10_baselined[ii,ichannel,:]  = (subject_data_ild10[ii,ichannel,:]  - np.nanmean(subject_data_ild10[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_ild10[:,:,0:int(5*fs)], axis=(0,1,2))
        
    #     subject_data_itd50_hbr_baselined[ii,ichannel,:]  = (subject_data_itd50_hbr[ii,ichannel,:]  - np.nanmean(subject_data_itd50_hbr[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_itd50_hbr[:,:,0:int(5*fs)], axis=(0,1,2))
    #     subject_data_itd500_hbr_baselined[ii,ichannel,:]  = (subject_data_itd500_hbr[ii,ichannel,:]  - np.nanmean(subject_data_itd500_hbr[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_itd500_hbr[:,:,0:int(5*fs)], axis=(0,1,2))
    #     subject_data_ild70n_hbr_baselined[ii,ichannel,:]  = (subject_data_ild70n_hbr[ii,ichannel,:]  - np.nanmean(subject_data_ild70n_hbr[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_ild70n_hbr[:,:,0:int(5*fs)], axis=(0,1,2))
    #     subject_data_ild10_hbr_baselined[ii,ichannel,:]  = (subject_data_ild10_hbr[ii,ichannel,:]  - np.nanmean(subject_data_ild10_hbr[:,ichannel,0:int(5*fs)], axis=(0,1)))#/np.nanstd(subject_data_ild10_hbr[:,:,0:int(5*fs)], axis=(0,1,2))


    #plot for this subject
    # fig, axes = plt.subplots(4, 2)
    # curr_ax = axes[0, 0]
    # curr_ax.plot(np.transpose(subject_data_itd50_baselined[ii, :, :]), 'r')
    # curr_ax = axes[1, 0]
    # curr_ax.plot(np.transpose(subject_data_itd50_baselined[ii, :, :]), 'r')
    # curr_ax = axes[2, 0]
    # curr_ax.plot(np.transpose(subject_data_itd50_baselined[ii, :, :]), 'r')
    # curr_ax = axes[3, 0]
    # curr_ax.plot(np.transpose(subject_data_itd50_baselined[ii, :, :]), 'r')
    
    # curr_ax = axes[0, 1]
    # curr_ax.plot(np.transpose(subject_data_itd50_hbr_baselined[ii, :, :]), 'b')
    # curr_ax = axes[1, 1]
    # curr_ax.plot(np.transpose(subject_data_itd50_hbr_baselined[ii, :, :]), 'b')
    # curr_ax = axes[2, 1]
    # curr_ax.plot(np.transpose(subject_data_itd50_hbr_baselined[ii, :, :]), 'b')
    # curr_ax = axes[3, 1]
    # curr_ax.plot(np.transpose(subject_data_itd50_hbr_baselined[ii, :, :]), 'b')
    
    # fig.suptitle(subject)
    # run a GLM to extract beta values for each condition

    # ---------------------------------------------------------------
    # -----------------     GLM                             ---------
    # ---------------------------------------------------------------

    # try to remove some of the conditions in the raw_haemo_filt annotations

    # raw_haemo_temp_crop = raw_haemo.crop(tmin=raw_haemo.annotations.onset[45] - 15)

    # raw_haemo_short_crop = get_short_channels(raw_haemo_temp_crop)
    # raw_haemo_filt_crop = get_long_channels(raw_haemo_temp_crop)

    # # try cropping
    # raw_haemo_filt_crop.resample(5)
    # raw_haemo_short_crop.resample(5)

    # raw_haemo_filt_crop.annotations.set_durations(glm_dur)

    # design_matrix_hbo = make_first_level_design_matrix(raw_haemo_filt_crop.pick(picks='hbo'),
    #                                                    drift_model=None,
    #                                                    high_pass=0.01,  # Must be specified per experiment
    #                                                    hrf_model='spm',
    #                                                    stim_dur=raw_haemo_filt_crop.annotations.duration)
    # # add_regs=filtered_signals)

    # design_matrix_hbo["Linear"] = np.arange(0, np.shape(design_matrix_hbo)[0])
    # design_matrix_hbo["ShortHbO"] = np.mean(raw_haemo_short_crop.copy().pick(picks="hbo").get_data(), axis=0)

    # design_matrix["ShortHbR"] = np.mean(short_chs.copy().pick(
    #     picks="hbr").get_data(), axis=0)
    # min_max_scaler = MinMaxScaler()
    # X_minmax = min_max_scaler.fit_transform(design_matrix_hbo)
    # design_matrix_min_max = pd.DataFrame(X_minmax, columns=design_matrix_hbo.columns.tolist())
    # if False:
    #     # plotting optional
    #     fig, ax1 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    #     ax_img = plot_design_matrix(design_matrix_min_max, ax=ax1)
    #     plt.show()
    #
    #     s = mne_nirs.experimental_design.create_boxcar(raw_haemo_filt, stim_dur=block_dur)
    #     plt.close()
    #     plt.plot(raw_haemo_filt.times[0:2000], s[0:2000, 3])
    #     plt.plot(design_matrix_hbo['Visual_Temporal'][0:1 + int(2000 / raw_haemo_filt.info['sfreq'])])
    #     plt.legend(["Stimulus", "Expected Response"])
    #     plt.xlabel("Time (s)")
    #     plt.ylabel("Amplitude")
    #     plt.show()

    # print(f'running GLM for subject {ii + 1}')

    # # pre-whiten
    # raw_haemo_filt_crop._data = np.subtract(raw_haemo_filt_crop._data,
    #                                         np.mean(raw_haemo_filt_crop._data, axis=1)[:, np.newaxis])
    # glm_est = run_glm(raw_haemo_filt_crop, design_matrix_hbo, noise_model='ar1')

    # record the glm est for each condition, for each subject
    # will adjust the beta values by the BH correction method

    # glm_est_df = glm_est.pick(picks='data', exclude='bads').to_dataframe()

    # # put into a larger array with all subjects data!
    # subject_data_itd50_GLM[ii, chan_indices_good] = glm_est_df.loc[glm_est_df['Condition'] == 'ild_0__itd_50']['theta']
    # subject_data_itd500_GLM[ii, chan_indices_good] = glm_est_df.loc[glm_est_df['Condition'] == 'ild_0__itd_500']['theta']
    # subject_data_ild70n_GLM[ii, chan_indices_good] = glm_est_df.loc[glm_est_df['Condition'] == 'ild_70n__itd_0']['theta']
    # subject_data_ild10_GLM[ii, chan_indices_good] = glm_est_df.loc[glm_est_df['Condition'] == 'ild_10__itd_0']['theta']



##############################
## Take mean during stim ####
############################

pfc_channels = []
stg_channels = []
# Take Means
index_stim_start = int(5*fs) # 8
index_stim_end = int(17.8*fs)
# ITD50
mean_during_stim_itd50 = np.nanmean(subject_data_itd50_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_itd50_hbr = np.nanmean(subject_data_itd50_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

# ITD500
mean_during_stim_itd500 = np.nanmean(subject_data_itd500_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_itd500_hbr = np.nanmean(subject_data_itd500_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

# ILD70n
mean_during_stim_ild70n = np.nanmean(subject_data_ild70n_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_ild70n_hbr = np.nanmean(subject_data_ild70n_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

#ILD10
mean_during_stim_ild10 = np.nanmean(subject_data_ild10_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_ild10_hbr = np.nanmean(subject_data_ild10_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

## Save breath uncorrected and corrected GLM data

# Uncorrected block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_hold.shape], names = names)
hold_df = pd.DataFrame({'subject_data_hold':subject_data_hold.flatten()},index=index)['subject_data_hold']
itd50_df = pd.DataFrame({'subject_data_itd50':subject_data_itd50_baselined.flatten()},index=index)['subject_data_itd50']
itd500_df = pd.DataFrame({'subject_data_itd500':subject_data_itd500_baselined.flatten()},index=index)['subject_data_itd500']
ild70n_df = pd.DataFrame({'subject_data_ild70n':subject_data_ild70n_baselined.flatten()},index=index)['subject_data_ild70n']
ild10_df = pd.DataFrame({'subject_data_ild10':subject_data_ild10_baselined.flatten()},index=index)['subject_data_ild10']
z = pd.concat([itd50_df,itd500_df,ild70n_df,ild10_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_uncorr_block_average__lateralization_target_{target_direction}_{masker_type}_masker.csv',index=True)

# Corrected block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_hold_bh_corr.shape], names = names)
hold_df = pd.DataFrame({'subject_data_hold_bh_corr':subject_data_hold_bh_corr.flatten()},index=index)['subject_data_hold_bh_corr']
itd50_df = pd.DataFrame({'subject_data_itd50_bh_corr':subject_data_itd50_bh_corr.flatten()},index=index)['subject_data_itd50_bh_corr']
itd500_df = pd.DataFrame({'subject_data_itd500_bh_corr':subject_data_itd500_bh_corr.flatten()},index=index)['subject_data_itd500_bh_corr']
ild70n_df = pd.DataFrame({'subject_data_ild70n_bh_corr':subject_data_ild70n_bh_corr.flatten()},index=index)['subject_data_ild70n_bh_corr']
ild10_df = pd.DataFrame({'subject_data_ild10_bh_corr':subject_data_ild10_bh_corr.flatten()},index=index)['subject_data_ild10_bh_corr']
z = pd.concat([itd50_df,itd500_df,ild70n_df,ild10_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_bh_corr_block_average_lateralization_target_{target_direction}_{masker_type}_masker.csv',index=True)


# Uncorrected block average means HBO
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd50.shape], names = names)
itd50_df = pd.DataFrame({'mean_during_stim_itd50':mean_during_stim_itd50.flatten()},index=index)['mean_during_stim_itd50']
itd500_df = pd.DataFrame({'mean_during_stim_itd500':mean_during_stim_itd500.flatten()},index=index)['mean_during_stim_itd500']
ild70n_df = pd.DataFrame({'mean_during_stim_ild70n':mean_during_stim_ild70n.flatten()},index=index)['mean_during_stim_ild70n']
ild10_df = pd.DataFrame({'mean_during_stim_ild10':mean_during_stim_ild10.flatten()},index=index)['mean_during_stim_ild10']
z = pd.concat([itd50_df,itd500_df,ild70n_df,ild10_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_lateralization_target_{target_direction}_{masker_type}_masker.csv',index=True)


# Uncorrected block average means HBR
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd50_hbr.shape], names = names)
itd50_df = pd.DataFrame({'mean_during_stim_itd50_hbr':mean_during_stim_itd50_hbr.flatten()},index=index)['mean_during_stim_itd50_hbr']
itd500_df = pd.DataFrame({'mean_during_stim_itd500_hbr':mean_during_stim_itd500_hbr.flatten()},index=index)['mean_during_stim_itd500_hbr']
ild70n_df = pd.DataFrame({'mean_during_stim_ild70n_hbr':mean_during_stim_ild70n_hbr.flatten()},index=index)['mean_during_stim_ild70n_hbr']
ild10_df = pd.DataFrame({'mean_during_stim_ild10_hbr':mean_during_stim_ild10_hbr.flatten()},index=index)['mean_during_stim_ild10_hbr']
z = pd.concat([itd50_df,itd500_df,ild70n_df,ild10_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_lateralization_target_{target_direction}_{masker_type}_masker_hbr.csv',index=True)


