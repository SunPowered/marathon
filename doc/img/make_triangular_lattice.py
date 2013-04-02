import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy

# we want the rotation matrices for a triangular lattice
precision = 3
def rotation_matrices():
	def Ry(alpha):
		cosa = numpy.around(numpy.cos(alpha), 2*precision)
		sina = numpy.around(numpy.sin(alpha), 2*precision)

		return numpy.array([[cosa, 0, -sina], [0, 1, 0], [sina, 0, cosa]])


	def Rz(alpha):
		cosa = numpy.around(numpy.cos(alpha),2*precision)
		sina = numpy.around(numpy.sin(alpha), 2*precision)

		return numpy.array([[cosa, -sina, 0], [sina, cosa, 0],[0, 0, 1]])

	step = numpy.pi/4
	matrices = [numpy.identity(3)]

	import ipdb; ipdb.set_trace()
	for ny in xrange(-2, 3):
		if ((ny == -2) or (ny == 2)):
			matrices.append(Ry(ny*step))
			continue
		for nz in xrange(8):
			if not ((nz == 0) and (ny == 0)):
				matrices.append(numpy.dot(Rz(nz*step), Ry(ny*step)))
	return matrices

matrices = rotation_matrices()

# Plot the unit branch

fig = plt.figure(1, figsize=[8,8])
ax = fig.add_subplot(111, projection="3d")

# Plot the flexible point
ax.scatter([0],[0],zs=[0], zdir="z", edgecolors='k', facecolor='yellow', s=500, marker='o', label="Flexible Point")

# Plot the branch atom
ax.scatter([1], [0], zs=[0], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

ax.plot([0,1], [0,0], [0,0], c='k', lw=2, alpha=1)

# Plot the remaining rotations
for R in matrices[1:]:
	rot_vec = numpy.dot(R, numpy.array([1,0,0]))
	ax.scatter(rot_vec[0], rot_vec[1], zs=rot_vec[2],
			edgecolors="k", facecolor="k", s=200, marker='o')
	ax.plot([0, rot_vec[0]], [0, rot_vec[1]], [0, rot_vec[2]], c="k", lw=2, alpha=0.5)

plt.title("Triangular Rotational Lattice")

limits = 1.1
ax.set_xlim([-limits, limits])
ax.set_ylim([-limits, limits])
ax.set_zlim([-limits, limits])
ax.plot([-limits, limits],[0,0], [0,0], ls='--', c='k' )
ax.plot([0,0], [-limits, limits], [0,0], ls='--', c='k' )
ax.plot([0,0], [0,0], [-limits, limits],ls='--', c='k' )
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#ax.legend()

# Try out the nifty XKCDify module
#import XKCDify

#ax = XKCDify.XKCDify(ax)
plt.show()
