import ScanImageTiffReader as scanimage
import os
import tifffile as tf
from skimage.transform import downscale_local_mean
import numpy as np
import pandas as pd

def extract_metadata(file_2p):
    filepath, filename = os.path.split(os.path.abspath(file_2p))
    
    all_meta = scanimage.ScanImageTiffReader(file_2p).metadata()
    all_metadata_file = filepath + os.path.sep + filename[:-4] + "_allMeta.txt"
    metadata_file = filepath + os.path.sep + filename[:-4] + "_Meta.txt"
    
    with open(all_metadata_file, 'w') as f:
        f.write(all_meta)
    
    all_meta = all_meta.splitlines() # split for easy extraction
    
    with open(metadata_file, 'w') as f:
        f.write(filepath + os.path.sep + filename + "\n\n")
        
        # Arm
        arm1_flag = [i for i in all_meta if 'Arm1' in i]
        arm2_flag = [i for i in all_meta if 'Arm2' in i]
        if len(arm1_flag) > 0 and len(arm2_flag) == 0:
            f.write('Arm 1\n\n')
        elif len(arm1_flag) == 0 and len(arm2_flag) > 0:
            f.write('Arm 2\n\n')
        else:
            f.write('Error in finding scan arm information from metadata\n\n')
        
        
        # Basic acquisition parameters
        f.write("Acquisition\n")
        laser_power = [i[i.find('[')+1:i.find('[')+3] for i in all_meta if 'powers ' in i][0]
        beamsplit = [i[i.find(']')-2:i.find(']')] for i in all_meta if 'powers ' in i][0]
        f.write("\tlaserPower = {}%\n".format(laser_power))
        f.write("\tbeamSplit = {}%\n".format(beamsplit))
        
        sample_position = [i[i.find('sample'):] for i in all_meta if 'samplePosition' in i][0]
        f.write("\t{}\n".format(sample_position))
        
        frame_rate = [i[i.find('scan'):] for i in all_meta if 'FrameRate' in i][0]
        f.write("\t{}\n\n".format(frame_rate))
        
        
        # Volumetric imaging parameters
        
        f.write("Volume\n")
        volume_rate = [i[i.find('scan'):] for i in all_meta if 'VolumeRate' in i][0]
        f.write("\t{}\n".format(volume_rate))
        
        num_slices = [i[i.find('numSlices'):] for i in all_meta if 'numSlices ' in i][0]
        f.write("\t{}\n".format(num_slices))
        
        step_size = [i[i.find('actual'):] for i in all_meta if 'actualStackZStepSize' in i][0]
        f.write("\t{}\n".format(step_size))
        
        start_pos = [i[i.find('stackZStartPos'):] for i in all_meta if 'stackZStartPos' in i][0]
        f.write("\t{}\n".format(start_pos))
        
        end_pos = [i[i.find('stackZEndPos'):] for i in all_meta if 'stackZEndPos' in i][0]
        f.write("\t{}\n\n".format(end_pos))
        
        
        # ROI information        
        roi_info = [i[i.find('['):i.find(']')+1] for i in all_meta if 'pixelResolutionXY' in i]
        if len(roi_info)>0:
            f.write("ROIs\n")
            for i, res in enumerate(roi_info):
                f.write("\troi{} PixelResolutionXY = {}".format(i, res) + "\n")
        
        return metadata_file
		

def segment_2p_video(file_2p, metadata_file, do_binning):    
    print("[+] Parsing " + metadata_file)
    num_rois = 0
    true_image_size = 0
    roi_y = []
        
    search_rois = True
    while search_rois:
        with open(metadata_file, 'r') as f:
            last_line = f.readlines()[-(num_rois+1)]
            if 'roi0 ' in last_line: # done searching for rois
                search_rois = False
        
        roi_y.insert(0, int(last_line[last_line.find(',')+1:-2]))
        true_image_size += roi_y[0]
        num_rois += 1
      
    print("\tFound {} ROIs".format(num_rois))
                    
    print("[+] Reading " + file_2p + "...")
    data2p = scanimage.ScanImageTiffReader(file_2p).data()
    num_junk_lines = data2p.shape[1] - true_image_size
    junk_lines_per_frame = int(num_junk_lines/num_rois)

    parent_directory = os.path.dirname(os.path.abspath(file_2p))

    for i, true_roi_height in enumerate(roi_y[::-1]): # index it backwards to account for junk lines
        roi_path = parent_directory + os.path.sep + 'roi_{}'.format(num_rois-i)
        set_path(roi_path)
        
        frame_height = true_roi_height + junk_lines_per_frame

        cropped_data = data2p[:, -frame_height:, :]
        data2p = data2p[:, 0:-frame_height, :]
		
        print("\t Writing region {} with no binning...".format(num_rois - i))
        tf.imwrite(roi_path + os.path.sep + 'data.tif', cropped_data, photometric='minisblack')

        if do_binning:
            print("\t\t[+] Binning region {}...".format(num_rois-i))
            set_path(roi_path)
            cropped_data = downscale_local_mean(cropped_data, (1, 2, 2)).astype('int16')
            print("\t\t[+] Writing region {} with binning...".format(num_rois - i))
            tf.imwrite(roi_path + os.path.sep + 'bin2x2x1' + os.path.sep + 'data.tif', cropped_data, photometric='minisblack')
		
        del cropped_data # memory management

        
		
def bin_2p(file2p, bin_factor=(1, 2, 2)):
    parent_directory = os.path.dirname(os.path.abspath(file2p))
    print("[+] Reading " + file2p + "...")
    data2p = scanimage.ScanImageTiffReader(file2p).data()
    data2p -= data2p.min()
    print("\t\t[+] Binning " + file2p + " by {}x{}x{}...".format(*bin_factor))
    binned_data = downscale_local_mean(data2p, bin_factor).astype('uint16')
    
    print("\t\t[+] Writing binned data...")
    set_path(parent_directory + os.path.sep + 'bin{}x{}x{}'.format(*bin_factor))
    tf.imwrite(parent_directory + os.path.sep + 'bin{}x{}x{}'.format(*bin_factor) + os.path.sep + 'data.tif', binned_data, photometric='minisblack')
    
def set_path(path):
    """
    Generates path recursively
    """
    if not os.path.exists(path):
        if len(path.split(os.path.sep)) == 1:
            raise Exception("Root path has been reached.")

        set_path(os.path.sep.join(path.split(os.path.sep)[:-1]))
        os.mkdir(path)