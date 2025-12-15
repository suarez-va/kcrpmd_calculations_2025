import sys
import os
import matplotlib.pyplot as plt   # plots
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np

from plot_utils import set_style

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

set_style()

inf = 10.
x_ar = np.linspace(-1.51, 1.51, 1000)
hnad1 = 1.
hnad2 = -1.
y_ar1 = (inf * (1 - np.heaviside(x_ar + 1.5, 0.5) + np.heaviside(x_ar - 1.5, 0.5))
         + hnad1 * (np.heaviside(x_ar + 0.5, 0.5) - np.heaviside(x_ar - 0.5, 0.5))) 
y_ar2 = (inf * (1 - np.heaviside(x_ar + 1.5, 0.5) + np.heaviside(x_ar - 1.5, 0.5))
         + hnad2 * (np.heaviside(x_ar + 0.5, 0.5) - np.heaviside(x_ar - 0.5, 0.5))) 

rp_0 = plt.imread("ring_polymers/rp_0.png")
rp_1 = plt.imread("ring_polymers/rp_1.png")
rp_kp1 = plt.imread("ring_polymers/rp_kp1.png")
rp_kp2 = plt.imread("ring_polymers/rp_kp2.png")
rp_kp3 = plt.imread("ring_polymers/rp_kp3.png")
rp_kp4 = plt.imread("ring_polymers/rp_kp4.png")
rp_kp5 = plt.imread("ring_polymers/rp_kp5.png")
rp_kp6 = plt.imread("ring_polymers/rp_kp6.png")
ib_0 = OffsetImage(rp_0, zoom=0.055)
ib_1 = OffsetImage(rp_1, zoom=0.055)
ib_kp1 = OffsetImage(rp_kp1, zoom=0.055)
ib_kp2 = OffsetImage(rp_kp2, zoom=0.055)
ib_kp3 = OffsetImage(rp_kp3, zoom=0.055)
ib_kp4 = OffsetImage(rp_kp4, zoom=0.055)
ib_kp5 = OffsetImage(rp_kp5, zoom=0.055)
ib_kp6 = OffsetImage(rp_kp6, zoom=0.055)
ab1_0 = AnnotationBbox(ib_0, (-1.0, 0.25), frameon=False)
ab1_1 = AnnotationBbox(ib_1, (1.0, 0.25), frameon=False)
ab1_kp1 = AnnotationBbox(ib_kp1, (0.0, 1.65), frameon=False)
ab1_kp2 = AnnotationBbox(ib_kp2, (-0.24, 1.25), frameon=False)
ab1_kp3 = AnnotationBbox(ib_kp3, (0.24, 1.25), frameon=False)
ab2_0 = AnnotationBbox(ib_0, (-1.0, 0.25), frameon=False)
ab2_1 = AnnotationBbox(ib_1, (1.0, 0.25), frameon=False)
ab2_kp4 = AnnotationBbox(ib_kp4, (0.0, -0.35), frameon=False)
ab2_kp5 = AnnotationBbox(ib_kp5, (-0.24, -0.75), frameon=False)
ab2_kp6 = AnnotationBbox(ib_kp6, (0.24, -0.75), frameon=False)

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(4.0, 4.0), sharex=True, sharey=False, dpi=300)
ax1.add_artist(ab1_0)
ax1.add_artist(ab1_1)
ax1.add_artist(ab1_kp1)
ax1.add_artist(ab1_kp2)
ax1.add_artist(ab1_kp3)
ax2.add_artist(ab2_0)
ax2.add_artist(ab2_1)
ax2.add_artist(ab2_kp4)
ax2.add_artist(ab2_kp5)
ax2.add_artist(ab2_kp6)
ax1.text(-0.1, 0.92, "(a)", transform=ax1.transAxes, fontsize=14)
ax2.text(-0.1, 0.92, "(b)", transform=ax2.transAxes, fontsize=14)
ax1.plot(x_ar, y_ar1, color='k', label='data1')
ax2.plot(x_ar, y_ar2, color='k', label='data2')
ax1.set_xlim(-1.6,1.6)
ax1.set_ylim(-0.2,2.0)
ax2.set_ylim(-1.3,0.9)
ax1.set_xticks([-1.5, -0.5, 0.5, 1.5])
ax1.set_yticks([])
ax2.set_yticks([])
ax2.set_xlabel("y")
ax1.set_ylabel(r"$F^{\mathrm{KC}}(y)$")
ax2.set_ylabel(r"$F^{\mathrm{KC}}(y)$")
ax1.set_title("")

plt.subplots_adjust(hspace=0.05, left=0.1, right=0.98, top=0.98, bottom=0.14)
plt.savefig('figure3.png')

