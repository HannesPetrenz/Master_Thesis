from scipy.io import savemat
import numpy as np
import glob
import os

npzFiles = glob.glob("L_track_barc.npz")
for f in npzFiles:
    fm = os.path.splitext(f)[0]+'.mat'
    d = np.load(f)
    savemat(fm, d)
    print('generated ', fm, 'from', f)