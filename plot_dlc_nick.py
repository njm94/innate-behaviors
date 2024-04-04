import os
import platform
import sys
from os.path import dirname, realpath

import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# from matplotlib import colors
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import itertools

filepath = realpath(__file__)
dir_of_file = dirname(filepath)
parent_dir_of_file = dirname(dir_of_file)
sys.path.append(parent_dir_of_file)
import tifffile as tif
import pandas as pd
import seaborn as sns
from pathlib import Path
import glob
from tqdm import tqdm
from scipy import signal
import scipy.io as sio
from scipy.ndimage import gaussian_filter
import imageio
import argparse
from configparser import ConfigParser
#import napari
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/numba_cache/'
import dlc2kinematics
# from pygifsicle import optimize
from get_list_data import expt1_data_list, expt2_data_list
from helper import *


# data_root = '/media/user/teamshare/nick/behavior/grooming/1p/'    
# data_list = expt1_data_list

data_root = '/media/user/teamshare/nick/behavior/grooming/2p/'    
data_list = expt2_data_list

# %%

for expt in tqdm(data_list):
    mouse_id, list_rec_dir = expt
    # if mouse id to process is passed, process only that mouse
    # out_dir = data_root + os.sep + mouse_id +  os.sep + "outputs"
    # set_path(out_dir)

    for rec_dir in tqdm(list_rec_dir):
        #if recording dir idc is passed, process only that dir
        try:
            dlc_file = glob.glob(data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + mouse_id + '_' + rec_dir + '*1030000.h5')[0]
        except:
            print('DLC tracking missing ', mouse_id + os.sep + rec_dir, ' continue to next')
            continue
        
        if not os.path.isfile(dlc_file.replace('.h5','.csv')):
            df_dlc, bodyparts, scorer = dlc2kinematics.load_data(dlc_file, smooth=True, filter_window=10, order=3)

            df_vel = dlc2kinematics.compute_velocity(df_dlc, bodyparts=['all'])
            df_acc = dlc2kinematics.compute_acceleration(df_dlc, bodyparts=['all'])
            df_speed = dlc2kinematics.compute_speed(df_dlc, bodyparts=['all'])   
            
            df_dlc.to_csv(dlc_file.replace('.h5','.csv'), index = False)
            df_vel.to_csv(dlc_file.replace('.h5','_vel.csv'), index = False)
            df_acc.to_csv(dlc_file.replace('.h5','_acc.csv'), index = False)
            df_speed.to_csv(dlc_file.replace('.h5','_speed.csv'), index = False)
        else:
            print('DLC csv files already created')

        