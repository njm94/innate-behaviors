#%%%% brain and behavior videos alignment based on blank frames

import sys
sys.path.append('../')
from os.path import dirname, realpath
filepath = realpath(__file__)
dir_of_file = dirname(filepath)
parent_dir_of_file = dirname(dir_of_file)
sys.path.append(parent_dir_of_file)
import tifffile as tif
import pandas as pd
from scipy import signal
from pathlib import Path
from tqdm import tqdm
from get_list_data import expt1_data_list, expt2_data_list, expt2_start_end_idx
from helper import *
import gc
import numpy as np
from scipy.io import savemat
from skimage.transform import rescale, resize, downscale_local_mean


# data_root = '/media/user/teamshare/nick/behavior/grooming/1p/'
data_root = '/media/user/teamshare/nick/behavior/grooming/2p/'

# data_root = '/mnt/njm/nick/behavior/grooming/1p/'
data_list = expt2_data_list
data_list_idx = expt2_start_end_idx

k = int(1000)

for ii, expt in enumerate(data_list):
    mouse_id, list_rec_dir = expt
    for jj, rec_dir in enumerate(tqdm(list_rec_dir)):
        # files = Path(data_root + os.sep + mouse_label + os.sep + d).glob('*.h264')
        beh_vid_file = data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + mouse_id + '_' + rec_dir + '.mp4'
        beh_log_file = data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + mouse_id + '_' + rec_dir + '.txt'
        if '1p' in data_root:
            brain_cam0_file = data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + 'cam0.tif'
            brain_frame_file = str(brain_cam0_file).replace('.tif', '_singleFrame.tif')
            brain_svd_file = str(brain_cam0_file).replace('.tif', '_svd.mat')
            brain_cam0_trim_file = str(brain_cam0_file).replace('.tif', '_trim.tif')


        beh_vid_trim_file = str(beh_vid_file).replace('.mp4', '_trim.mp4')
        beh_log_trim_file = str(beh_log_file).replace('.txt', '_trim.txt')
        
        
        beh_svd_file = str(beh_vid_file).replace('.mp4', '_svd.mat')
        beh_motion_file = str(beh_vid_file).replace('.mp4', '_ME.avi')
        beh_motion_svd_file = str(beh_vid_file).replace('.mp4', '_MEsvd.mat')

        if '1p' in data_root and os.path.isfile(brain_svd_file):
            print(mouse_id + "_" + rec_dir + " already processed. Skipping...")
            continue
        
        # if os.path.isfile(beh_vid_trim_file) and \
        #     (not os.path.isfile(brain_cam0_file) or os.path.isfile(brain_svd_file)) and \
        #         (os.path.isfile(beh_svd_file)) and \
        #             os.path.isfile(beh_motion_svd_file):
        #     # print('Session already aligned, Skip ',beh_vid_file)
        #     continue

        if True:#not os.path.isfile(beh_vid_trim_file):
            if '1p' in data_root:
                beh_start, beh_stop = get_dark_frames(beh_stack)
            elif '2p' in data_root: # issues with TTL to scanimage cause inconsistent event timings: use manually selected indices
                beh_start, beh_stop = data_list_idx[ii][jj]
                print(beh_start, beh_stop)

            
            print('Processing ', beh_vid_file)
            beh_stack =  open_cv_read_video(beh_vid_file)
            beh_log = pd.read_csv(beh_log_file, sep='\t', parse_dates=['time'])

            
                
            beh_log_trim = beh_log.iloc[beh_start:beh_stop]
            beh_log_trim.frame = beh_log_trim.frame - beh_start
            beh_log_trim.to_csv(beh_log_trim_file, index=False, sep='\t')
            beh_stack_trim = beh_stack[beh_start:beh_stop]
            open_cv_write_video(beh_vid_file, beh_vid_trim_file, beh_start, beh_stop)
        
            # memory management
            del beh_log
            del beh_log_trim
            del beh_stack
        else:
            print("Reading trim behavior video")
            beh_stack_trim = open_cv_read_video(beh_vid_trim_file)
        
        num_beh_frames = beh_stack_trim.shape[0]

        print("Downsampling frames")
        beh_stack_trim = np.array(beh_stack_trim[::5,:,:,0])

        print("Doing local mean")
        beh_stack_trim = downscale_local_mean(beh_stack_trim, (1, 3, 3)).astype('uint8')

        if not os.path.isfile(beh_svd_file):
            beh_stack_for_svd = np.reshape(beh_stack_trim, newshape=(beh_stack_trim.shape[0], int(beh_stack_trim.shape[1]*int(beh_stack_trim.shape[2])))).T
            print("Behavior video SVD")
            U, s, V = movie_svd(beh_stack_for_svd, k)
            V_upsampled = signal.resample(V, num_beh_frames, axis=1)
            mdic = {"U": U, "s": s, "V": V_upsampled}
            savemat(beh_svd_file, mdic)
            del U
            del V_upsampled
            del V
            del s
            del beh_stack_for_svd
        
        if not os.path.isfile(beh_motion_svd_file) or not os.path.isfile(beh_motion_file):
            beh_stack_ME = calculate_ME(beh_stack_trim)
            del beh_stack_trim
            open_cv_write_video_from_arr(beh_motion_file, beh_stack_ME, fps=90)

            beh_stack_ME = np.reshape(beh_stack_ME, newshape=(beh_stack_ME.shape[0], int(beh_stack_ME.shape[1]*int(beh_stack_ME.shape[2])))).T
            print("Behavior video SVD")
            U, s, V = movie_svd(beh_stack_ME, k)
            V_upsampled = signal.resample(V, num_beh_frames, axis=1)
            mdic = {"U": U, "s": s, "V": V_upsampled}
            savemat(beh_motion_svd_file, mdic)
            del U
            del V_upsampled
            del V
            del s
            del beh_stack_ME
            
        #    U, V = movie_svd(beh_stack_ME, k)
        #    V_upsampled = signal.resample(V, num_beh_frames, axis=1)
        #    mdic = {"U": U, "s": s, "V": V_upsampled}
        #    savemat(beh_motion_svd_file, mdic)

        gc.collect()        

        if not os.path.isfile(brain_cam0_file):
            print('No Brain data, Skip ', brain_cam0_file)
            continue
        
        # if os.path.isfile(brain_cam0_trim_file):
        #     print('Trim file already computed, Skip ', brain_cam0_trim_file)
        #     continue

        if not os.path.isfile(brain_svd_file) or not os.path.isfile(brain_frame_file):
            print("reading the brain stack")
            brain_stack = tif.imread(brain_cam0_file)
            brain_start, brain_stop = get_dark_frames(brain_stack)
            brain_stack_trim = brain_stack[brain_start:brain_stop]
            del brain_stack
            tif.imwrite(brain_frame_file, brain_stack_trim[0, :, :], photometric='minisblack')

            brain_stack_trim = np.reshape(brain_stack_trim, newshape=(brain_stack_trim.shape[0], int(brain_stack_trim.shape[1] * brain_stack_trim.shape[2]))).T
            print("SVD on brain")
            U, s, V = movie_svd(brain_stack_trim, k)
            V_upsampled = signal.resample(V, num_beh_frames, axis=1)            
            mdic = {"U": U, "s": s, "V": V_upsampled}
            savemat(brain_svd_file, mdic)

        # brain_stack_trim = signal.resample(brain_stack_trim, num_beh_frames, axis=0).astype(np.uint16)
        # tif.imwrite(brain_cam0_trim_file, brain_stack_trim, photometric='minisblack')



# %%
