import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Plot the unit branch

fig = plt.figure(1, figsize=[8,8])
ax = fig.add_subplot(111, projection="3d")

# Plot the flexible point
ax.scatter([0],[0],zs=[0], zdir="z", edgecolors='k', facecolor='yellow', s=500, marker='o', label="Flexible Point")

# Plot the branch atoms
ax.scatter([1], [0], zs=[0], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

ax.scatter([0], [0], zs=[-1], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

ax.scatter([0], [1], zs=[0], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

ax.scatter([-1], [0], zs=[0], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

ax.scatter([0], [-1], zs=[0], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

ax.scatter([0], [0], zs=[1], zdir="z", edgecolors="k", facecolor="k", s=200, marker='o', label="Branch Atom")

# Plot the bonds

ax.plot([0,1], [0,0], [0,0], c='k', lw=2, alpha=1)
ax.plot([0,0], [0,0], [0,-1], c='k', lw=2, alpha=0.5)
ax.plot([0,0], [0,1], [0,0], c='k', lw=2, alpha=0.5)
ax.plot([0,-1], [0,0], [0,0], c='k', lw=2, alpha=0.5)
ax.plot([0,0], [0,-1], [0,0], c='k', lw=2, alpha=0.5)
ax.plot([0,0], [0,0], [0,1], c='k', lw=2, alpha=0.5)




plt.title("Cubic Rotational Lattice")

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
