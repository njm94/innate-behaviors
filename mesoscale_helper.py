import numpy as np
import cv2
import matplotlib.pyplot as plt
import os






def calculate_df_f0(frames):
    """
    Calculate df/f0, the fractional change in intensity for each pixel
    and the variance of df/f0
    """
    frames = frames.astype(np.float32)
    baseline = np.mean(frames, axis=0)
    df_f0 = np.divide(
        np.subtract(frames, baseline), baseline
    )

    return df_f0


def calculate_df_f0_moving(frames, n=144, axis=0):
    """
    Calculate df/f0, the fractional change in intensity for each pixel
    and the variance of df/f0 with a moving baseline (f0)
    """
    frames = frames.astype(np.float32)
    baseline = np.cumsum(frames, axis, dtype=np.float32)
    baseline[n:] = baseline[n:] - baseline[:-n]
    baseline[n-1:] = baseline[n-1:] / n

    prelim = np.arange(n)

    # df_f0[:n-1] = df_f0[:n-1] / prelim[1:]
    baseline[:n - 1] = baseline[:n - 1] / prelim[1:, None, None]

    # df_f0 = frames - baseline

    df_f0 = np.divide(
        np.subtract(frames, baseline), baseline
    )
    # del frames, baseline
    # df_f0[
    #     numpy.where(numpy.isnan(df_f0))
    # ] = -1  # Make the nans black.

    return df_f0



class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        self.slices, rows, cols = X.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.X[self.ind, :, :])
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[self.ind, :, :])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()



def load_frames(filename, color):
    """
    Load frames of .h264/5 as color channel(s) or B&W frames as numpy array

    :param filename: path to video file
    :type: str
    :param color: one of ('red', 'green', 'blue', 'all', False), with False for B&W
    :type: str or bool

    :return: video frames
    :type: numpy.ndarray
    """
    channel = {"red": 0, "green": 1, "blue": 2, False:'', "all":''}.get(color)
    if channel is None:
        raise AttributeError(
            "Argument 'color' must be one of ('red', 'green', 'blue', 'all', False)"
        )
    cap = cv2.VideoCapture(filename)
    num_frames = 0
    first_frame = True
    while cap.isOpened:
        ret, frame = cap.read()
        if not ret:
            break
        if first_frame:
            height, width, _ = frame.shape
            first_frames = False
        num_frames +=1
    cap.release()
    if num_frames == 0:
        raise Exception(f"No frames found in 'filename'")
    if color == 'all':
        frames = np.zeros((num_frames, height, width, 3), dtype=np.uint8)
    else:
        frames = np.zeros((num_frames, height, width), dtype=np.uint8)
    index = 0
    cap = cv2.VideoCapture(filename)
    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break
        if not color:
            frames[index] = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        elif color == 'all':
            frames[index] = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        else:
            frames[index] = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)[
                ..., channel
            ]
        index+=1
    cap.release()
    return frames


def set_path(path):
    """
    Generates path recursively
    """
    if not os.path.exists(path):
        if len(path.split(os.path.sep)) == 1:
            raise Exception("Root path has been reached.")

        set_path(os.path.sep.join(path.split(os.path.sep)[:-1]))
        os.mkdir(path)