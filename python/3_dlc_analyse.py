#%%%% DeepLabCut analyse the behavior .mp4 videos
import os
from tqdm import tqdm
import glob
from get_list_data import expt1_data_list, expt2_data_list
import deeplabcut as dlc

# config_path = '/home/pankaj/teamshare/TM_Lab/nick/behavior/grooming/mouse_grooming-nick-2023-11-27/config.yaml'
config_path = '/media/user/teamshare/nick/behavior/grooming/mouse_grooming-nick-2023-11-27/config.yaml'

# data_root = '/home/pankaj/teamshare/TM_Lab/nick/behavior/grooming/1p/'
data_root = '/media/user/teamshare/nick/behavior/grooming/2p/'
# data_root = '/media/user/teamshare/nick/behavior/grooming/1p/'

# data_list = expt1_data_list
data_list = expt2_data_list
# data_list.extend(expt2_data_list)

### To export the model
# dlc.export_model(config_path, iteration=None, shuffle=1, trainingsetindex=0, snapshotindex=None, TFGPUinference=True, overwrite=False, make_tar=False)

for expt in data_list:
    mouse_id, list_rec_dir = expt

    for rec_dir in tqdm(list_rec_dir):
        beh_vid_file = data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + mouse_id + '_' + rec_dir + '_trim.mp4'
        print(beh_vid_file)
        if os.path.isfile(beh_vid_file):
            print(glob.glob(data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + mouse_id + '_' + rec_dir + '*1030000.h5'))
            print("\n\n")
            if len(glob.glob(data_root + os.sep + mouse_id + os.sep + rec_dir + os.sep + mouse_id + '_' + rec_dir + '*1030000.h5')) == 0:
                print("hello")
                
                dlc.analyze_videos(config_path, [beh_vid_file], videotype='.mp4', shuffle=1, trainingsetindex=0, save_as_csv = False)
                # dlc.plot_trajectories(config_path, [beh_vid_file], videotype='.mp4', showfigures=False)
                # dlc.create_labeled_video(config_path, beh_vid_file, videotype='.mp4', draw_skeleton=False, overwrite=True)

                # dlc.analyzeskeleton(config_path, beh_vid_file, videotype='.mp4', shuffle=1, trainingsetindex=0, save_as_csv=True)
