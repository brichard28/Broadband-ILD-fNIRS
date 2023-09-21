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
all_fnirs_data_folders = ['C:/Users/benri/Nextcloud/nirs/data/2023-07-20/2023-07-20_002',
'C:/Users/benri/Nextcloud/nirs/data/2023-07-20/2023-07-20_003',
'C:/Users/benri/Nextcloud/nirs/data/2023-07-26/2023-07-26_001',
'C:/Users/benri/Nextcloud/nirs/data/2023-08-02/2023-08-02_001']



subject_ID = ['NDARYZ656HJ9','NDARCD778KPR','NDARMY829TKN','NDARLU426TBZ']
curr_subject_ID = ['NDARYZ656HJ9','NDARCD778KPR','NDARMY829TKN','NDARLU426TBZ']


def individual_analysis(fnirs_data_folder, ID):
    
    os.environ["OMP_NUM_THREADS"] = "1"

    raw_intensity = read_raw_nirx(fnirs_data_folder, verbose=False)
    
    raw_intensity.annotations.rename({'15.0':'Masker Noise ITD 50 Targ Left',
                                      '14.0':'Masker Noise ILD 10 Targ Right',
                                      '13.0':'Masker Noise ILD 10 Targ Left',
                                      '12.0':'Masker Speech ILD 10 Targ Right',
                                      '11.0':'Masker Speech ITD 500 Targ Left',
                                      '10.0':'Masker Noise ITD 500 Targ Left',
                                      '9.0':'Masker Speech ITD 50 Targ Left',
                                      '8.0':'Masker Noise ITD 500 Targ Right',
                                      '7.0':'Masker Speech ITD 50 Targ Right',
                                      '6.0':'Masker Speech ILD 10 Targ Left',
                                      '5.0':'Masker Speech ITD 500 Targ Right',
                                      '4.0':'Masker Noise ITD 50 Targ Right',
                                      '3.0':'Hold',
                                      '2.0':'Exhale',
                                      '1.0':'Inhale'})
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
    
    # View Consistency across channels
    fig2, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 6))
    clims = dict(hbo=[-10, 10], hbr=[-10, 10])
    epochs["Masker Speech ITD 500 Targ Left"].average().plot_image(axes=axes[:, 0], clim=clims)
    epochs["Masker Speech ITD 500 Targ Right"].average().plot_image(axes=axes[:, 1], clim=clims)
    for column, condition in enumerate(["Masker Speech ITD 500 Targ Left", "Masker Speech ITD 500 Targ Right"]):
        for ax in axes[:, column]:
            ax.set_title("{}: {}".format(condition, ax.get_title()))
            
    # Plot standard evoked
    evoked_dict = {
    "left/HbO": epochs["Masker Speech ITD 500 Targ Left"].average(picks="hbo"),
    "left/HbR": epochs["Masker Speech ITD 500 Targ Left"].average(picks="hbr"),
    "right/HbO": epochs["Masker Speech ITD 500 Targ Right"].average(picks="hbo"),
    "right/HbR": epochs["Masker Speech ITD 500 Targ Right"].average(picks="hbr"),
    }

    # Rename channels until the encoding of frequency in ch_name is fixed
    for condition in evoked_dict:
        evoked_dict[condition].rename_channels(lambda x: x[:-4])

    color_dict = dict(HbO="#AA3377", HbR="b")
    styles_dict = dict(right=dict(linestyle="dashed"))
           
  

    
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
    mne.viz.plot_evoked_topo(
        epochs["Masker Speech ITD 500 Targ Left"].average(picks="hbo"), color="b", axes=axes, legend=True
        )
    mne.viz.plot_evoked_topo(
        epochs["Masker Speech ITD 500 Targ Right"].average(picks="hbo"), color="r", axes=axes, legend=True
        )
    
    # Save Block Averages
    epochs_this_subject = epochs.to_data_frame()
    epochs_this_subject.to_csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\" + sub + " block averages.csv")
    
    return raw_haemo, epochs



    
    

for sub in curr_subject_ID:  

    # Create path to file based on experiment info
    which_folder = subject_ID.index(sub)
    fnirs_data_folder = all_fnirs_data_folders[which_folder]

    # Analyse data and return both ROI and channel results
    raw_haemo, epochs = individual_analysis(fnirs_data_folder, sub)



