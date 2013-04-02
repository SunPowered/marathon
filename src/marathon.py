#!/usr/bin/env python
"""
marathon.py

A python script to rotate molecules parsed in  Protein Database files (.pdb) 
format.  A molecule is parsed, the flexible points are calculated, and these
serve as the rotation joints.  90 degree and 45 degree rotations are 
available, which rotate the bond at each loop by 90/45 degree increments in 3d, 
excepting the incoming bond direction.  


The script is run from a command prompt through the python interpreter

$>python marathon -o output_dir molecule.pdb 


All options are given by calling the program's help option.

$>python marathon.py --help


The output directory is set with the -o/--output options, expecting a following
argument of a location to save new .pdb files and optionally rendered strucure
plots.

Cubic rotations (i.e. 90 deg iterations) are set by default.  Triangular
rotations (i.e. 45 degree iterations) are set with the -t flag.

The script is capable of plotting the structure of the molecule and the 
rotations using the -p/--plot flag.  The plots are saved by default in .pdf
format in a subdirectory of the output director called "plots".

The -i/--interactive flag indicates that all rotation plots should be displayed
on the screen in an interactive mode, allowing structure exploration and atom
inspection using a mouse


PDB file support gracefully borrored from mmLib <http://pymmlib.sourceforge.net/>


--------
TODO
-------
Read some options from settings file
"""

import sys
import os
import copy
import shutil

from itertools import product
from pprint import pprint

try:
	import numpy
except ImportError:
	print "This program requires numpy to be installed.  Try `pip install numpy`"
	sys.exit(1)

# For now, use the module provided, but maybe in future merge this module in here
import PDB

__version__ = "0.9.1"

###############################################################
# Program Options
##############################################################
precision = 3 	# Number of decimals to save in float values
plot_extension = ".pdf" # Default plot extension to use
figure_size = [10, 10] #Figure size in inches.  Maybe put this in the matplotlibrc file?
sep="-"
rmsd_fname = "rmsd.txt" # The default name to use for the file saving the RMSD of a molecule
##############################################################

class BondInterferenceError(Exception):
	pass

def Rz(alpha):
	"""Rotation matrix around the z axis"""

	cosa = numpy.around(numpy.cos(alpha),2*precision)
	sina = numpy.around(numpy.sin(alpha), 2*precision)

	return numpy.array([[cosa, -sina, 0], [sina, cosa, 0],[0, 0, 1]])

def Ry(alpha):
	"""Rotation matrix around the y axis"""
	cosa = numpy.around(numpy.cos(alpha), 2*precision)
	sina = numpy.around(numpy.sin(alpha), 2*precision)

	return numpy.array([[cosa, 0, -sina], [0, 1, 0], [sina, 0, cosa]])


def rotation_matrices(N=4, verbose=False):
	"""Split a unit circle into N equal sections and permute rotatations around
	the Y and Z axis this many times.  I.e. cubic rotations, N=4.  Triangular 
	rotations, N=8.  
	
	Assumes the original bond vector is along the X axis, these rotations must
	be adjusted if the original vector is off
	"""
	if N not in [4, 8]:
		raise ValueError("Only N of 4 and 8 are currently supported")

	step = 2*numpy.pi/N
	
	if verbose:
		print "Rotation Matrices with N: {}".format(N)
		print "step: {}".format(step)

	
	# Its nice to have the identity  matrix come first
	matrices = [numpy.identity(3)]

	
	for ny in xrange(-N/4, N/4+1):
		if (ny == -N/4) or (ny == N/4):
			if verbose:
				print "Ny: {}".format(ny)
			matrices.append(Ry(ny*step))
			continue
		for nz in xrange(N):
			if (nz == 0) and (ny==0):
				continue
			'''
			if (nz*step == numpy.pi) and (ny == 0 ):
				# Ignore the rotation that turns it on the incoming bond
				if verbose:
					print "Ignoring nz {} / ny {}".format(nz, ny)
				continue
			'''
			if verbose:
				print "Nz: {} / Ny: {}".format(nz,ny)

			matrices.append(numpy.dot(Rz(nz*step), Ry(ny*step)))
	return matrices
	
class PDBMolecule(object):
	"""PDB molecule data structure.  This holds all information 
	pertaining to one molecule.  It is capable of reading from a 
	.pdb file, scanning the molecule for flexible points and 
	parsing branches from these loops, rotating branches
	as well as writing the rotated structure back to a file.
	
	Can also plot the molecule interactively or to file"""

	def __init__(self, file_path, verbose=False):
		"""Arguments: 
			file_name:  The file to read
			verbose: 	A boolean to print more detailed info to the console"""
		self.file_path= file_path
		self.file_name = os.path.split(self.file_path)[1]
		self.verbose = verbose
		self.data = []
		self.atoms = []
		self.flexible_points = []
		self.branches = []
		self.rotated_branches = []
		self.read()

	def __str__(self):
		return str(self.data)
	
	def plot_molecule(self, file=None, title=None):
		"""Plot the current molecular structure.
		
		The coordinates of each atom are placed with
		approrpriate colouring and size along with
		all bonds.  
		
		If the file argument is provided than save 
		plot to a file, otherwise plot interactively.
		
		Optionally set the title of the plot, defaults
		to the original .pdb file name"""

		try:
			import matplotlib.pyplot as plt
			from mpl_toolkits.mplot3d import Axes3D
		except ImportError:
			print "Cannot import matplotlib"
			return

		def atom_onpick(event):
			"""Selection event for atoms.  
			
			currently just display the atom info to the console"""
			lbl = event.artist.get_label()
			ind = event.ind
			
			print
			print "Atom Selection"
			pprint(plot_atoms[lbl]["atoms"][ind[0]])

	
		plot_title = title or self.file_name
		
		# get axes
		fig = plt.figure(1, figsize=figure_size)
		ax = fig.add_subplot(111, projection="3d")

		# Extract the various atoms of interest as well as all the others
		carbon_atoms = [atom for atom in self.atoms if atom.element=="C"]
		flexible_points = [atom for atom in self.flexible_points]
		nitrogen_atoms = [atom for atom in self.atoms if atom.element == "N" and atom not in flexible_points] 
		other_atoms = [atom for atom in self.atoms if (atom not in carbon_atoms) and (atom not in nitrogen_atoms) and (atom not in flexible_points)]


		# atom plotting styles
		default_style = {"color": "k",
				"size": 50,
				"edgecolors": 'none'}
		carbon_style = default_style.copy()
		nitrogen_style = default_style.copy()
		nitrogen_style.update({"color": 'r'})
		il_style = default_style.copy()
		il_style.update({"color": "w",
			"size": 100,
			"edgecolors": "k"})
		other_style = default_style.copy()
		other_style.update({"color": 'y'})
	
		# structure to store the bond start and end points
		lines = {"x":[], "y":[], "z":[] }

		# Build the bond lines from the atoms
		for atom in self.atoms:
			for bond_id in atom.bonds:
				b_atom = self.get_atom_by_id(bond_id)
				lines["x"].append([atom.X, b_atom.X])
				lines["y"].append([atom.Y, b_atom.Y])
				lines["z"].append([atom.Z, b_atom.Z])

		
		leg = {"handlers":[], "labels":[]}
		plot_atoms = {
			"Carbon": {"atoms": carbon_atoms,
				"style" : carbon_style},
			"Nitrogen": {"atoms": nitrogen_atoms,
				"style": nitrogen_style},
			"Internal Loop": {"atoms": flexible_points,
				"style": il_style},
			"Other": {"atoms": other_atoms, 
				"style": other_style}}

		# add each group of atoms to the plot
		for label, itm in plot_atoms.iteritems():

			xs = [atom.X for atom in itm["atoms"]]
			ys = [atom.Y for atom in itm["atoms"]]
			zs = [atom.Z for atom in itm["atoms"]]

			ax.scatter(xs, ys, zs=zs, 
					zdir="z", edgecolors=itm["style"]['edgecolors'], 
					facecolor=itm["style"]["color"], s=itm["style"]["size"], 
					marker="o", label=label,
					picker=True)

			leg["handlers"].append(plt.Circle((0,0), fc=itm["style"]["color"]))
			leg["labels"].append(label)
		
		# display the legend for the atoms
		ax.legend(leg["handlers"], leg["labels"], title="Atoms")

		# connect selection event with custom handler
		fig.canvas.mpl_connect("pick_event", atom_onpick)

		# add the bond lines
		for xs, ys, zs in zip(lines["x"], lines["y"], lines["z"]):
			ax.plot(xs, ys, zs, c="black", lw=1)
	
		# set the plot axes labels and titles
		ax.set_xlabel("X")
		ax.set_ylabel("Y")
		ax.set_zlabel("Z")
		plt.title("{}".format(plot_title))

		# save or show?
		if file is not None:
			plt.savefig(file, bbox=0.0)
		else:
			plt.show()
		plt.clf()
		plt.close()

	def load_atoms(self):
		"""Take atoms read from the PDB module and extract relevant info for our needs"""
		if not self.data:
			if self.verbose:
				print "No data found for the molecule"
			return

		# Quick and dirty way of extracting the atoms 
		# and bonds in the molecule
		for atom in self.data:
			if "element" in atom:
				# add atom
				if self.verbose:
					print "Found atom: {}".format(atom)
				self.atoms.append(PDBAtom(atom))
			elif "serialBond1" in atom:
				# Add bond
				if self.verbose:
					print "Found bond: {}".format(atom)
				self.add_bond(atom)
	
	def get_flexible_points(self):
		"""Flexible points defined as any element N, surrounded by C's
		
		Saves the flexible points to the molecule, does not return anything"""
		joints = []
		if self.verbose:
			print "Looking for flexible points"

		for atom in self.atoms:
			if atom.element != "N":
				continue

			bonds = atom.bonds
			if len(bonds) < 2:
				continue

			is_IL = True
			for bond in bonds:
				for a in self.atoms:
					if a.seq_id == bond:
						if a.element != "C":
							is_IL = False
							break

			if is_IL:
				self.flexible_points.append(atom)
				[self.branches.append([atom, bond]) for bond in atom.bonds]

		if self.verbose:
			print "Found {} flexible points and {} branches".format(len(self.flexible_points), len(self.branches))

	def rotate_branch(self, branch_id, R):
		"""Rotates a branch by a specified rotation matrix. The 
		branch is given as a joint_atom and a bond_atom pair.
		The rotation matrix is applied to all atoms around the 
		original branch bond direction, this is always relative
		to the bond before rotation"""
		
		
		# Get the branch
		branch = self.branches[branch_id]
		
		joint_atom = branch[0]
		branch_atom = self.get_atom_by_id(branch[1])
		
		# normalize rotation_axis
		rotation_axis = branch_atom.coordinates - joint_atom.coordinates
		rotation_axis = rotation_axis / numpy.sqrt(numpy.dot(rotation_axis, rotation_axis))

		# This is the direction the bond axis will point after rotation 
		rotated_axis = numpy.around(numpy.dot(R, rotation_axis), precision)

		# Ensure that this rotated bond axis is not parallel to another bond axis at 
		# The flexible joint
		for bond in joint_atom.bonds:
			# Ignore itself
			if bond == branch_atom.seq_id:
				continue
			if bond not in self.rotated_branches:
				continue
			bond_axis = self.get_atom_by_id(bond).coordinates - joint_atom.coordinates
			bond_axis = bond_axis / numpy.sqrt(numpy.dot(bond_axis, bond_axis))
			bond_axis = numpy.around(bond_axis, precision)
			if numpy.allclose(numpy.dot(bond_axis, rotated_axis), 1):
				# Parallel vectors
				raise BondInterferenceError("Rotation of branch {} will collide with bond {}".format(branch_id, bond))

		# If no collision, continue with the rotation
		# Get the subtree
		st = self.get_subtree(joint_atom, branch_atom)
		
		# Rotate all atoms of all sub branches
		for atom in st:
			atom.rotate(R, joint_atom.coordinates, rotation_axis)
		self.rotated_branches.append(branch_atom.seq_id)
		
	def center_coordinates(self):
		"""Translate the structure so that an atom in the middle lies at 0,0"""
		
		# Get the average x,y,z values for the structure
		avgs = [
				numpy.mean([atom.X for atom in self.atoms]),
				numpy.mean([atom.Y for atom in self.atoms]),
				numpy.mean([atom.Z for atom in self.atoms])]

		# find the closest distanced atom to this point
		min_i = 999999
		min_dist = 999999
		for i, atom in enumerate(self.atoms):
			distance_to_avg = numpy.abs(atom.coordinates - avgs)
			if distance_to_avg < min_dist:
				min_i = i
				min_dist = distance_to_avg

		center_atom = self.atoms[min_i].coordinates

		# center the structure around this point
		for atom in self.atoms:
			atom.X = atom.X - center_atom[0]
			atom.Y = atom.Y - center_atom[1]
			atom.Z = atom.Z - center_atom[2]

	def update_data(self):
		"""Update the coordinates of the data elements"""
		for atom in self.atoms:
			for a in self.data:
				if not "element" in a:
					# Not an atom item
					continue
				if a['resSeq'] == atom.seq_id:
					a['x'] = atom.X
					a['y'] = atom.Y
					a['z'] = atom.Z

	def get_atom_by_id(self, atom_id):
		"""Search the available atoms to find one named by its id"""
		r =  [atom for atom in self.atoms if atom.seq_id == atom_id]
		if len(r) == 0:
			return None
		else:
			return r[0]

	def get_subtree(self, atom, child):
		"""Returns the molecule branch after a given atom, traversing towards a child"""
		data=[]	
		if len(child.bonds) == 1:
			return [child]

		for bond_id in child.bonds:
			if bond_id == atom.seq_id:
				continue

			bonded_atom = self.get_atom_by_id(bond_id)
			sub = [child]
			sub.extend(self.get_subtree(child, bonded_atom))
			data.extend(sub)
			
		return data


	def add_bond(self, bond):
		"""Add a bond connection"""
		n1 = bond.get("serial")
		n2 = bond.get("serialBond1")
		[atom.bonds.append(n2) for atom in self.atoms if atom.seq_id == n1]
		[atom.bonds.append(n1) for atom in self.atoms if atom.seq_id == n2]

	def read(self):
		"""Read a .pdb file"""
		if not os.path.isfile(self.file_path):
			raise ValueError("File {} not found".format(self.file_path))
		self.data = PDB.PDBFile()
		self.data.load_file(open(self.file_path, 'r'))
		self.load_atoms()
		self.get_flexible_points()
		if self.verbose:
			print "Loaded file {}".format(self.file_path)
	

	def write(self, file_name):
		"""write current molecule to a .pdb file"""
		self.update_data()
		if self.verbose:
			print "Writing to file {}".format(file_name)

		self.data.save_file(open(file_name, 'w'))

	@property
	def rmsd(self):
		"""Calculates RMSD of molecule.
		
		The square root of the mean of the square of the distance between 
		all atoms i, j where i != j in the molecule"""
		return  numpy.sqrt(numpy.mean([(atom.coordinates - atom_p.coordinates)**2 for atom in self.atoms for atom_p in self.atoms if atom_p != atom ]))


class PDBAtom(object):
	"""Data structure for each atom, holds the atom informationn as well as 
	bonding information"""

	def __init__(self, atom):
		"""Atom is a dict read from the PDB file"""
		self.element = atom.get("element")
		self.seq_id = atom.get("resSeq", 0)
		self.X = atom.get("x", 0.0)
		self.Y = atom.get("y", 0.0)
		self.Z = atom.get("z", 0.0)
		

		self.bonds = []

	@property
	def coordinates(self):
		return numpy.array([self.X, self.Y, self.Z])

	def __repr__(self):
		return "PDBAtom {} ({}) at ({}, {}, {}).  Bond({})".format(self.seq_id, self.element, self.X, self.Y, self.Z, self.bonds)

	def rotate(self, matrix, point=numpy.array([0,0,0]), direction=numpy.array([1, 0,0]), verbose=False):
		"""Rotates the  atom with respect to a point and direction
		matrix is a 3x3 rotation matrix
		point is a 3 dimensional vector
		direction is a 3 dimensional unit vector in the base coordinates
		"""
		if verbose:
			print "Rotating atom with rotation matrix"
			pprint(matrix)
			print "rotation point: {}".format(point)
			print "direction: {}".format(direction)

		# Normalize direction
		direction = direction / numpy.sqrt(numpy.dot(direction, direction))
		rc = self.coordinates - point
		rc_mag = numpy.dot(rc, rc)
		if numpy.allclose(rc_mag, 0):
			if verbose:
				print "rotation around self detected"
				return
				
		# Normalize 
		rc_norm = rc/numpy.sqrt(rc_mag)
		


		# Coordinates must be rotated to be in line with the parent
		# This rotation is along the vector perpendicular to both rc and direction
		rot_unit = numpy.cross(rc_norm, direction)
		dot = numpy.dot(rc_norm, direction)
		
		rot_theta = numpy.arccos(numpy.dot(rc_norm, direction))
		rot_R = rotation_matrix(direction, rot_theta)
		rc = numpy.dot(rot_R, rc)

		
		# apply rotation matrix to this rotated point
		rc = numpy.dot(matrix, rc)

		# Revert back to original coordinates
		rc = numpy.dot(numpy.transpose(rot_R), rc)
		rc = rc + point

		# Save coordinates to atom
		self.X = round(rc[0], precision)
		self.Y = round(rc[1], precision)
		self.Z = round(rc[2], precision)


def rotation_matrix(axis,theta):
	"""Euler-Rodrigues angles"""
	axis = axis/numpy.sqrt(numpy.dot(axis,axis))

	a = numpy.cos(theta/2)
	b,c,d = -axis*numpy.sin(theta/2)
	return numpy.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
		

def rotation_permutations_from_file(pdb_file, rotation="cubic", verbose=False, 
		plot=False, out_dir=None, interactive=False, rmsd=False, detailed=False):

	"""The main function to process files.
	
	This takes a single pdb file and performs the necessary rotations on each
	flexible point.  All combinations are saved to a file and optionally plotted
	
	The rotation type is specified with the 'rotation' argument.
	verbose=True spits out more information to the console
	plot=True plots each rotation iteration
	out_dir specify the output directory to save files
	interactive=True Dont save plots but direct them to the display"""

	print "Running rotational permutations on file:  {}".format(pdb_file)
	
	# Read file and find flexible points

	pdb_orig = PDBMolecule(pdb_file, verbose=verbose)

	base_name = os.path.splitext(pdb_orig.file_name)[0]

	# Check the output directory
	if not out_dir:
		out_dir = os.path.dirname(pdb_file)
	file_out_dir = os.path.join(out_dir, base_name)

	if not os.path.isdir(file_out_dir):
		if verbose:
			print "Output Directory not found.  Creating: {}".format(file_out_dir)
		os.makedirs(file_out_dir)
	
	plot_out_dir = os.path.join(file_out_dir, "plots")
	if not os.path.isdir(plot_out_dir) and plot:
		os.mkdir(plot_out_dir)

	# Delete the rmsd file if it already exists
	rmsd_filename = os.path.join(file_out_dir, rmsd_fname)
	if rmsd:
		if os.path.isfile(rmsd_filename):
			if verbose:
				print "Removing pre existing {} file".format(rmsd_fname)
			os.remove(rmsd_filename)
		rmsd_fd = open(rmsd_filename, 'w')
		rmsd_fd.write("Iteration\tRMSD\n")

	# Get the rotations for this round
	if rotation == "triangular":
		Rs = rotation_matrices(N=8)
	else:
		Rs = rotation_matrices(N=4)


	# For all unique branches 
	# Compute permutations with the R rotation matrices at loop
	branches = pdb_orig.branches
	if verbose:
		print "A total of {} points and {} branches to permute around".format(len(pdb_orig.flexible_points), len(branches))
	
	#For each branch, get the rotation permutations 
	perms_storage = []
	rot_perms = product(xrange(len(Rs)), repeat=len(branches))
		
	# Iterate over each set of rotations
	perm_count = 1	
	for rot_perm in rot_perms:

		f = copy.deepcopy(pdb_orig)
		try:
			if detailed:
				rot_name = sep.join(["B{}R{}".format(N, R) for N, R  in enumerate(rot_perm)])
			else:
				rot_name = str(perm_count)
			
			file_name = rot_name + sep + base_name

			# Rotate all the branches with their respective rotation matrix for this iteration
			for N, R in enumerate(rot_perm):
				f.rotate_branch(N, Rs[R])

		except BondInterferenceError:
			if verbose:
				print "Bond Interference detected.  Ignoring rotation"
		else:
			# calculate rmsd?
			if rmsd:
				rmsd_fd.write("{}\t{}\n".format(rot_name, f.rmsd))

			# write the file and plot it
			f.write(os.path.join(file_out_dir, file_name+".pdb"))
			if plot:
				if interactive:
					plot_file = None
				else:
					plot_file = os.path.join(plot_out_dir, file_name + plot_extension)
				f.plot_molecule(file=plot_file, title=file_name)
			perm_count += 1
		finally:
			del f

	if verbose:
		print "Completed permutations on file {}".format(pdb_file)	
	if rmsd:
		rmsd_fd.close()

def main():
	"""Runs the script parsing arguments"""
	
	import argparse
	parser = argparse.ArgumentParser(description="Program to parse a PDB file, identify isolation loops, and permute molecular rotations around those loops and write back to a set of PDB files")
	parser.add_argument("args", help="One or more filenames or directories", nargs="+")
	parser.add_argument("-v", "--verbose", help="Print details to console", action="store_true")
	parser.add_argument("-o", "--output", help="Output new PDB files to this directory", default=os.getcwd(), action="store")
	parser.add_argument("-c", "--cubic", help="Rotate around a cubic structure (default), i.e. 90 deg", action="store_true")
	parser.add_argument("-t", "--triangular", help="Roatate around a triangular structure, i.e. 45 deg", action="store_true")
	parser.add_argument("-p", "--plot", help="Plot the rotated molecules in a `plots` subfolder", default=False, action="store_true")
	parser.add_argument("-i", "--interactive", help="Plot figures interactively", action="store_true")
	parser.add_argument("-r", "--rmsd", help="Calculate the root means square distance of each iteration and save all values to a file", action="store_true", default=False)
	parser.add_argument('-d', '--detailed', help="Save the iteration names with detailed information for each branch and rotation number, otherwise just use the iteration counter as a name", default=False, action="store_true")
	args = parser.parse_args()

	print
	print "Permuting Rotations v{}".format(__version__)
	print 
	
		
	if args.verbose:
		print "Verbose enabled"

		#print "Arguments received: {}".format( args)

	print "Saving output files to directory {}".format(args.output)	
	if not os.path.isdir(args.output):
		if args.verbose:
			print "Creating {}".format(args.output)
		os.makedirs(args.output)

	rotation_method = "cubic"
	if args.triangular:
		rotation_method = "triangular"
	if args.verbose:
		print "Using {} rotation structure".format(rotation_method)

	if args.rmsd:
		if args.verbose:
			print "Calulating RMSD for each iteration and saving results to a sub  directory of {}".format(args.output)

	if args.interactive:
		args.plot = True
	
	# Check the positional arguments
	if args.args:
		for arg in args.args:
			if os.path.isfile(arg):
			
				rotation_permutations_from_file(arg, 
						rotation=rotation_method, 
						verbose=args.verbose, 
						plot=args.plot, 
						out_dir=args.output, 
						interactive=args.interactive,
						rmsd=args.rmsd,
						detailed=args.detailed)

			elif os.path.isdir(arg):
				# Its a directory, permute each file individually
				for f in os.listdir(arg):
			
					# Run permutations on files in the directory
					pdb_file = 	os.path.join(arg, f)
					rotation_permutations_from_file(pdb_file, 
							rotation=rotation_method, 
							verbose=args.verbose, 
							plot=args.plot, 
							out_dir=args.output, 
							interactive=args.interactive, 
							rmsd=args.rmsd,
							detailed=args.detailed)

			else:
				print "Unknown argument: {}.  Skipping.".format(arg)
	print "All Done.  Have a nice day"


if __name__ == "__main__":
	main()	
