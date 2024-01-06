import ScanImageTiffReader as scanimage
import os
import tifffile as tf
from skimage.transform import downscale_local_mean
import numpy as np
import pandas as pd
import cv2
import numpy as np
import os
import ffmpeg
import numpy as np
from scipy import stats
from tqdm import tqdm


def calculate_ME(arr):
    print("[+] Calculating motion energy")
    for i, frame in enumerate(tqdm(arr)):
        if i > 0:
            arr_ME[i-1] = cv2.absdiff(frame, old_frame)
        else:
            arr_ME = np.empty(shape=(arr.shape[0]-1, arr.shape[1], arr.shape[2]), dtype=np.uint8)
        old_frame = frame
    
    return arr_ME


def extract_metadata(file_2p):
    print("[+] Extracting metadata from 2p video")
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
		

def ffmpeg_read_video(beh_vid_file):
    probe = ffmpeg.probe(beh_vid_file)
    video_info = next(s for s in probe['streams'] if s['codec_type'] == 'video')
    width = int(video_info['width'])
    height = int(video_info['height'])
    fps = eval(video_info['avg_frame_rate'])
    # num_frames = int(video_info['nb_frames'])
    out, err = (
        ffmpeg
            .input(beh_vid_file)
            .output('pipe:', format='rawvideo', pix_fmt='rgb24', loglevel="quiet")
            .run(capture_stdout=True)
    )
    beh_stack = (
        np
            .frombuffer(out, np.uint8)
            .reshape([-1, height, width, 3])
    )
    return beh_stack, fps


def ffmpeg_video_info(beh_vid_file):
    probe = ffmpeg.probe(beh_vid_file)
    video_info = next(s for s in probe['streams'] if s['codec_type'] == 'video')
    width = int(video_info['width'])
    height = int(video_info['height'])
    fps = eval(video_info['avg_frame_rate'])
    # num_frames = int(video_info['nb_frames'])

    return width, height, fps


def ffmpeg_write_video(rgb_stack, video_path, fps=30):
    process = (
        ffmpeg
            .input('pipe:', format='rawvideo', pix_fmt='rgb24', framerate=fps, s='{}x{}'.format(rgb_stack.shape[2], rgb_stack.shape[1]))
            .output(video_path, pix_fmt='yuv420p')
            .overwrite_output()
            .run_async(pipe_stdin=True)
    )
    process.stdin.write(
        rgb_stack
            .tobytes()
    )
    process.stdin.close()
    process.wait()

     
def get_dark_frames(frames):
    start_index = 0
    end_index = len(frames)

    mean_zs = stats.zscore(frames.mean(axis=tuple(range(1, frames.ndim)))) #mean along all axis except first
    ixs = np.where(mean_zs[:100]>-5)[0]
    start_index = ixs[0] if len(ixs) > 0 else 0
    ixs = np.where(mean_zs[-100:]<-5)[0]
    end_index = (len(mean_zs) -100 + ixs[0]) if len(ixs) > 0 else len(mean_zs)

    return (start_index, end_index)


def get_dark_frames_gradient_method(behaviour_frames, sigma=15, show_plot=False, spacetime=False):
    if spacetime:
        means = np.mean(behaviour_frames, axis=1)
    else:
        means = np.mean(np.mean(behaviour_frames, axis=1), axis=1)
    grad = np.gradient(means)
    mean = np.mean(grad)
    std = np.std(grad)

    if show_plot:
        plt.figure()
        plt.plot(np.gradient(means))
        plt.show()

    indeces = grad[np.where(np.abs(grad) > mean + std * sigma)]
    start = np.where(grad == indeces[0])
    stop = np.where(grad == indeces[-1])
    assert (len(start) == 1)
    assert (len(stop) == 1)
    return (start[0][0], stop[0][0])
    
    
def movie_svd(movie, k):
    print("[+] Computing SVD")
    U, s, Vt = np.linalg.svd(movie, full_matrices=False)
    U_k = U[:, :k]
    s_k = np.diag(s[:k])
    Vt_k = Vt[:k, :]

    return U_k, s_k, Vt_k


def open_cv_read_video(video_path):
    cap = cv2.VideoCapture(video_path)
    frameCount = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    frameWidth = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    frameHeight = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    buf = np.empty((frameCount, frameHeight, frameWidth, 3), np.dtype('uint8'))
    for i in tqdm(range(frameCount)):
        ret, buf[i] = cap.read()
        if not ret:
            break

    cap.release()
    return buf

    
def open_cv_write_video(video_path, output_path, beh_start, beh_stop):
    cap = cv2.VideoCapture(video_path)
    n_frames_old = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    fps = int(cap.get(cv2.CAP_PROP_FPS))
    
    fourcc = cv2.VideoWriter_fourcc(*'MP4V')
    out = cv2.VideoWriter(output_path, fourcc, fps, (width,height))               

    counter = 0
    for i in tqdm(range(n_frames_old)):
        ret, frame = cap.read()
        if not ret:
            break

        if counter >= beh_start and counter < beh_stop:
            out.write(frame)
        counter += 1
        
    cap.release()
    out.release()


def open_cv_write_video_from_arr(output_path, arr, fps=30):
    print("[+] Writing video")
    height, width = arr[0].shape[:2]
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    out = cv2.VideoWriter(output_path, fourcc, fps, (width, height), isColor=False)
    for frame in tqdm(arr):
        out.write(frame)
        
    out.release()


def reconstruct_svd(U, s, V):
    print("[+] Reconstructing movie from SVD")
    return np.dot(U, np.dot(s, V))

   
def segment_2p_video(file_2p, metadata_file, do_binning):    
    print("[+] Parsing ROIs from metadata")
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
                    
    print("[+] Reading 2-photon data")
    data2p = scanimage.ScanImageTiffReader(file_2p).data()
    print(data2p.dtype)
    
    if data2p.max() < np.iinfo(np.int16).max:
        data2p -= data2p.min()
        data2p = data2p.astype('uint16')
    
    num_junk_lines = data2p.shape[1] - true_image_size
    junk_lines_per_frame = int(num_junk_lines/num_rois)

    parent_directory = os.path.dirname(os.path.abspath(file_2p))

    for i, true_roi_height in enumerate(roi_y[::-1]): # index it backwards to account for junk lines
        roi_path = parent_directory + os.path.sep + 'roi_{}'.format(num_rois-i)
        set_path(roi_path)
        
        frame_height = true_roi_height + junk_lines_per_frame

        cropped_data = data2p[:, -frame_height:, :]
        data2p = data2p[:, 0:-frame_height, :]


        if do_binning:
            print("[+] Binning region {}".format(num_rois-i))
            roi_path = roi_path + os.path.sep + 'bin2x2x1'
            set_path(roi_path)
            cropped_data = downscale_local_mean(cropped_data, (1, 2, 2)).astype('uint16')

        roi_filename = roi_path + os.path.sep + 'data.tif'
        print("[+] Writing region {}".format(num_rois - i))
        tf.imwrite(roi_filename, cropped_data, photometric='minisblack')


def set_path(path):
    """
    Generates path recursively
    """
    if not os.path.exists(path):
        if len(path.split(os.path.sep)) == 1:
            raise Exception("Root path has been reached.")

        set_path(os.path.sep.join(path.split(os.path.sep)[:-1]))
        os.mkdir(path)


def xt(data, fs):
    return np.divide(range(0, data.shape[0]), fs)
