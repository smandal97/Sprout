import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import h5py

def make_sphere(center, arrsize, radius):
    distance = np.linalg.norm(np.subtract(np.indices(arrsize).T, center), axis=-1)
    mask = np.ones(arrsize) * (distance < radius)
    return mask

f = h5py.File('/Users/mandal0/Desktop/Numerical/Sprout/checkpoint_0006.h5', 'r')

i0f = 0.5; j0f = 0.5; k0f = 0.;
i0 = int(i0f*f['Data']['Cells'][0,...].shape[2]) 
j0 = int(j0f*f['Data']['Cells'][0,...].shape[1])
k0 = int(k0f*f['Data']['Cells'][0,...].shape[0])

data = f['Data']['Cells'][5,...]

rho = f['Data']['Cells'][0,k0,...]
p   = f['Data']['Cells'][1,k0,...]
vx  = f['Data']['Cells'][2,k0,...]
vy  = f['Data']['Cells'][3,k0,...]
vz  = f['Data']['Cells'][4,k0,...]
ps  = f['Data']['Cells'][5,k0,...]
T   = f['Grid']['T'][0]


dx = f['Grid']['dxs'][0]   ; dy = f['Grid']['dxs'][1]   ; dz = f['Grid']['dxs'][2]
Ox = f['Grid']['Origin'][0]; Oy = f['Grid']['Origin'][1]; Oz = f['Grid']['Origin'][2]
xl = Ox; yl = Oy; zl = Oz
xr = Ox + f['Data']['Cells'][0,...].shape[2]*dx
yr = Oy + f['Data']['Cells'][0,...].shape[1]*dy
zr = Oy + f['Data']['Cells'][0,...].shape[0]*dz

plt.figure(figsize=[6,5])
plt.title("T=%.3e" %T)
plt.imshow(np.log10(rho), origin='lower', extent=[xl,xr,yl,yr], cmap='inferno_r')
plt.colorbar(fraction=0.1)
plt.show()
