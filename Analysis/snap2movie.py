import os
import numpy as np
import matplotlib.pyplot as plt

#parameters
n = 400
r_w = 1e-2
ip_dir = '/Users/mandal0/Desktop/Numerical/Sprout_tools/kb_snaps/'
op_dir = '/Users/mandal0/Desktop/Numerical/Sprout_tools/kb_snaps/'
n_dpi = 200

#make images
for i in range(n):
    fname  = 'slicemap_'+f'{i:04d}'+".dat"
    img    = np.loadtxt(ip_dir+fname)
    header = np.loadtxt(ip_dir+fname, dtype=str, max_rows=1, comments=None)
    dl     = float(header[1])
    t      = float(header[2])
    x_w    = np.linspace(0.,img.shape[0]*dl,1000)
    y_w    = np.sqrt(r_w**2.-x_w**2.)
    with plt.style.context('dark_background'):
        plt.xlim([0,img.shape[0]*dl])
        plt.ylim([0,img.shape[0]*dl])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.imshow(np.log10(img), extent=[0,img.shape[0]*dl,0,img.shape[0]*dl], origin='lower', cmap='inferno')
        plt.plot(x_w,y_w,'r--')
        plt.title("$\\mathrm{log}_{10}\\rho$ (t=%.2e)" %t)
        #plt.show()
        plt.savefig(op_dir+'im_'+f'{i:04d}'+'.jpg', dpi=n_dpi)


#make movie
shell_cmd = "cd " + op_dir + "; pwd; ffmpeg -framerate 12 -pattern_type glob -i '*.jpg' video.mp4"
stream = os.popen(shell_cmd)
output = stream.read()
output