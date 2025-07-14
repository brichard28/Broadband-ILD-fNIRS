# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:06:19 2024

@author: benri
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
                          data_root + '2024-07-01/2024-07-01_001',
                          data_root + '2024-07-01/2024-07-01_002',
                          data_root + '2024-07-02/2024-07-02_001',
                          data_root + '2024-07-02/2024-07-02_002']


data = read_raw_nirx(f"{all_fnirs_data_folders[0]}/", verbose=False, preload=True)

# ---------------------------------------------------------------
# -----------------      Visualize Montage           ---------
# ---------------------------------------------------------------  
subjects_dir = str(mne.datasets.sample.data_path()) + "/subjects"
#mne.datasets.fetch_hcp_mmp_parcellation(subjects_dir=subjects_dir, accept=True)
labels = mne.read_labels_from_annot(
    "fsaverage", "HCPMMP1", "lh", subjects_dir=subjects_dir
)
labels_combined = mne.read_labels_from_annot(
    "fsaverage", "HCPMMP1_combined", "lh", subjects_dir=subjects_dir
)
view_map = {
"left-lat": np.r_[17, 18, 20, 21],
"left-frontal": np.r_[1,2,4,5 ],
"right-frontal": np.r_[7,9],
"right-lat": np.r_[11, 12, 14, 15],
}

ch_names_dict = {'S1_D1': '1',
 'S1_D2': '2',
 'S1_D8': '',
 'S2_D1': '3',
 'S2_D2': '4',
 'S2_D9': '',
 'S3_D3': '5',
 'S3_D10': '',
 'S4_D3': '6',
 'S4_D11': '',
 'S5_D4': '7',
 'S5_D5': '8',
 'S5_D12': '',
 'S6_D4': '9',
 'S6_D5': '10',
 'S6_D13': '',
 'S7_D6': '11',
 'S7_D7': '12',
 'S7_D14': '',
 'S8_D6': '13',
 'S8_D7': '14',
 'S8_D15': ''}

#fig_montage = mne_nirs.visualisation.plot_3d_montage(
#    data.info, view_map=view_map, subjects_dir=subjects_dir, src_det_names=defaultdict(lambda: ''), ch_names=ch_names_dict)

#fig_montage.set_tight_layout(tight=True)
picks = mne.pick_types(data.info, meg=False, fnirs=True)
dists = mne.preprocessing.nirs.source_detector_distances(data.info,picks=picks)
raw_filtered = data.pick(picks[dists >= 0.01])
info_filtered = raw_filtered.info
brain = mne.viz.Brain("fsaverage", subjects_dir=subjects_dir, background="w", cortex="0.5")
brain.add_sensors(info_filtered, trans="fsaverage", fnirs=["channels","pairs"])