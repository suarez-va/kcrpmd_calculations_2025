import sys
import os
import h5py
import matplotlib.pyplot as plt   # plots
import numpy as np

from liblibra_core import *
from libra_py import units

from plot_utils import set_style, add_hbar, add_abar

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

set_style()


#_sys_3_method_3_fix_s_leps_-7.60e-03_hw_1
#_sys_3_method_3_fix_y_leps_-7.60e-03_hw_-1

with h5py.File("../_sys_3_method_3_fix_s_leps_-7.60e-03_hw_1/libra_data/_reptraj_pt1/mem_data.hdf", 'r') as f:
    R_fix_s_pt1 = f["q/data"][:, 0, :]
    y_fix_s_pt1 = f["y_aux_var/data"][:, 0]

with h5py.File("../_sys_3_method_3_fix_s_leps_-7.60e-03_hw_1/libra_data/_reptraj_pt2/mem_data.hdf", 'r') as f:
    R_fix_s_pt2 = f["q/data"][:, 0, :]
    y_fix_s_pt2 = f["y_aux_var/data"][:, 0]

with h5py.File("../_sys_3_method_3_fix_y_leps_-7.60e-03_hw_-1/libra_data/_reptraj_pt1/mem_data.hdf", 'r') as f:
    R_fix_y_pt1 = f["q/data"][:, 0, :]
    y_fix_y_pt1 = f["y_aux_var/data"][:, 0]

with h5py.File("../_sys_3_method_3_fix_y_leps_-7.60e-03_hw_-1/libra_data/_reptraj_pt2/mem_data.hdf", 'r') as f:
    R_fix_y_pt2 = f["q/data"][:, 0, :]
    y_fix_y_pt2 = f["y_aux_var/data"][:, 0]

R_fix_s_pt2 = np.flip(R_fix_s_pt2, axis=0)
R_fix_y_pt2 = np.flip(R_fix_y_pt2, axis=0)

R_fix_s = np.concatenate((R_fix_s_pt2, R_fix_s_pt1), axis=0)
R_fix_y = np.concatenate((R_fix_y_pt2, R_fix_y_pt1), axis=0)

s_fix_s = R_fix_s[:,0]
q_fix_s = R_fix_s[:,-1]
s_fix_y = R_fix_y[:,0]
q_fix_y = R_fix_y[:,-1]
print("-------")
print(s_fix_s)
print("-------")
print(q_fix_s)
print("-------")
print(s_fix_y)
print("-------")
print(q_fix_y)


exit()
print("-------")
print(R_fix_y.shape)
exit()
print(np.flip(R_fix_y_pt2, axis=1))
print(R_fix_y_pt1.shape)
print(R_fix_y_pt2.shape)


exit()


fig, ax = plt.subplots(dpi=200)
#ax.set_xlim(-2.1,1.1)
#ax.set_ylim(-21.0,-9.0)
#ax.set_xticks([-2, -1, 0, 1])
#ax.set_yticks([-20, -18, -16, -14, -12, -10])
#ax.set_xlabel(r"log(β$K_0$)", fontsize = 15)
#ax.set_ylabel(r"log($k_{\mathrm{ET}}$)", fontsize = 15)
ax.set_title("")
#ax.legend(ncol=2, fontsize=10, loc='upper left', handletextpad=0.2,  columnspacing=0.9)

#plt.subplots_adjust(left=0.18, right=0.98, top=0.98, bottom=0.17)
plt.savefig('fig7.png')

