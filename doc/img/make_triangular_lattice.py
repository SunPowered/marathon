"""
This module constructs a diagram displaying the cubic rotations of a vector

There must be a 3D projection axis, There must be clear vectors coming from the 
origin along the 90 degree increments.  The vectors must be labeled """

from matplotlib import rc
font_styles = {'family':'sans-serif',
		'sans-serif':'Verdana'
		}
rc('font', **font_styles)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'src'))

try:
	import marathon
except ImportError:
	print "canot find marathon"
	exit()

node_styles = dict(edgecolors='k', 
					facecolor='w', 
					s=300, 
					marker='o')
line_styles = dict(lw=3, c='k')
other_styles = dict(text_offset=1.5)

triangular_vecs = marathon.triangular_bond_directions()

def plot_bond_vector(ax, vector, label=""):
	ax.scatter(vector[0], vector[1], vector[2], **node_styles)
	ax.plot([0, vector[0]], [0, vector[1]], [0, vector[2]], **line_styles)
	offset = other_styles.get('text_offset', 2.0)
	#t_vec = offset * vector
	ax.text(offset*vector[0], offset*vector[1], offset*vector[2], label, 
			fontdict=dict(weight='bold', size=14))

def plot_unit_cube(ax):
	cube_vecs = [
			[[-1, 1], [1, 1], [1,1]], 
			[[1, 1], [-1, 1], [1,1]], 
			[[1, 1], [1, 1], [-1,1]], 
			[[-1, 1], [-1, -1], [1,1]], 
			[[1, 1], [-1, 1], [-1,-1]], 
			[[-1, -1], [1, 1], [-1,1]], 
			[[-1, 1], [1, 1], [-1,-1]], 
			[[-1, -1], [-1, 1], [1,1]], 
			[[1, 1], [-1, -1], [-1,1]],
			[[-1, 1], [-1, -1], [-1,-1]], 
			[[-1, -1], [-1, 1], [-1,-1]], 
			[[-1, -1], [-1, -1], [-1,1]], 
			] 
	for x, y, z in cube_vecs:
		ax.plot(x, y, z, c='k', ls='--', lw=3, alpha=0.7)


fig = plt.figure(1, figsize=[8,8])
ax3d = fig.add_subplot(111, projection="3d")

for N, vec in enumerate(triangular_vecs):
	rot_label = "R{}".format(N)
	plot_bond_vector(ax3d, vec, label=rot_label)


plot_unit_cube(ax3d)
plt.title("Triangular Rotational Lattice", fontdict=dict(weight='bold', size=18))

limits = 2.5
axis_lbl_styles=dict(style='italic', size=10)
axis_lbl_offset = 1.1
ax3d.set_xlim([-limits, limits])
ax3d.set_ylim([-limits, limits])
ax3d.set_zlim([-limits, limits])
ax3d.plot([-limits, limits],[0,0], [0,0], ls='--', c='k' )
ax3d.text(limits*axis_lbl_offset , 0, 0, 'X', fontdict=axis_lbl_styles)
ax3d.plot([0,0], [-limits, limits], [0,0], ls='--', c='k' )
ax3d.text(0, limits*axis_lbl_offset , 0, "Y", fontdict=axis_lbl_styles)
ax3d.plot([0,0], [0,0], [-limits, limits],ls='--', c='k' )
ax3d.text(0, 0, limits*axis_lbl_offset,  "Z", fontdict=axis_lbl_styles)
ax3d.set_axis_off()

# Try out the nifty XKCDify module
#import XKCDify

#ax = XKCDify.XKCDify(ax)
plt.show()
