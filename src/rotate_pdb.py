#!/usr/bin/env python
"""
PDB support files gracefully borrored from mmLib <http://pymmlib.sourceforge.net/>
"""
import sys
import os
import copy
from pprint import pprint

try:
	import numpy
except ImportError:
	print "This program requires numpy to be installed.  Try `pip install numpy`"
	sys.exit(1)

import PDB
#import rotations

__version__ = 0.5
precision = 2 	# Number of decimals to save in float values
plot_extension = ".pdf"
figure_size = [10, 10]

def Rz(alpha):
	"""Rotation matrix around the z axis"""

	cosa = numpy.cos(alpha)
	sina = numpy.sin(alpha)

	return numpy.array([[cosa, -sina, 0], [sina, cosa, 0],[0, 0, 1]])

def Ry(alpha):
	"""Rotation matrix around the y axis"""
	cosa = numpy.cos(alpha)
	sina = numpy.sin(alpha)

	return numpy.array([[cosa, 0, -sina], [0, 1, 0], [sina, 0, cosa]])


def rotation_matrices(N=4, verbose=False):
	"""Split the circle into N equal sections and permute rotatations around
	the Y and Z axis this many times.  I.e. cubic rotations, N=4.  Triangular 
	rotations, N=8.  Ignore the case when the rotation angle is [pi/2, 0]
	
	Assumes the original bond vector is along the X axis, these rotations must
	be adjusted if the original vector is off
	"""
	if N not in [4, 8]:
		raise ValueError("Only N of 4 and 8 are currently supported")

	step = 2*numpy.pi/N
	
	if verbose:
		print "Rotation Matrices with N: {}".format(N)
		print "step: {}".format(step)

	
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
			
			if (nz*step == numpy.pi) and (ny == 0 ):
				# Ignore the rotation that turns it on the incoming bond
				if verbose:
					print "Ignoring nz {} / ny {}".format(nz, ny)
				continue
			
			if verbose:
				print "Nz: {} / Ny: {}".format(nz,ny)


			matrices.append(numpy.dot(Rz(nz*step), Ry(ny*step)))
	return matrices
	
class PDBfile(object):
	"""PDB File to read modify and write"""

	def __init__(self, file_path, verbose=False):
		"""Arguments: 
			file_name:  The file to read"""
		self.file_path= file_path
		self.file_name = os.path.split(self.file_path)[1]
		self.verbose = verbose
		self.data = []
		self.atoms = []
		self.internal_loops = []
		self.branches = []
		self.read()

	def __str__(self):
		return str(self.data)
	
	def plot_molecule(self, file=None):
		try:
			import matplotlib.pyplot as plt
			from mpl_toolkits.mplot3d import Axes3D
		except ImportError:
			print "Cannot import matplotlib"
			return

		# get axes
		fig = plt.figure(1, figsize=figure_size)
		ax = fig.add_subplot(111, projection="3d")


		# plot carbons
		carbon_atoms = [atom for atom in self.atoms if atom.element=="C"]
		internal_loops = [atom for atom in self.internal_loops]
		nitrogen_atoms = [atom for atom in self.atoms if atom.element == "N" and atom not in internal_loops] 
		other_atoms = [atom for atom in self.atoms if (atom not in carbon_atoms) and (atom not in nitrogen_atoms) and (atom not in internal_loops)]

		lines = {"x":[], "y":[], "z":[] }

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
	
		
		for atom in self.atoms:
			for bond_id in atom.bonds:
				b_atom = self.get_atom_by_id(bond_id)
				lines["x"].append([atom.X, b_atom.X])
				lines["y"].append([atom.Y, b_atom.Y])
				lines["z"].append([atom.Z, b_atom.Z])

		
		#ax.plot(lines["x"], lines["y"], zs=lines["z"], lw=2, c="black")
		#ax.plot(x,y,z, lw=2, c="black")

		leg = {"handlers":[], "labels":[]}
		plot_atoms = {
			"Carbon": {"atoms": carbon_atoms,
				"style" : carbon_style},
			"Nitrogen": {"atoms": nitrogen_atoms,
				"style": nitrogen_style},
			"Internal Loop": {"atoms": internal_loops,
				"style": il_style},
			"Other": {"atoms": other_atoms, 
				"style": other_style}}

		def atom_onpick(event):
			lbl = event.artist.get_label()
			ind = event.ind
			
			print
			print "Selection"
			#print "lbl: {} / ind: {}".format(lbl, ind)
			pprint(plot_atoms[lbl]["atoms"][ind[0]])

			#atom = plot_atoms[lbl]["atoms"][ind[0]]
			#print "Selected atom: {}".format(atom)

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
			
		ax.legend(leg["handlers"], leg["labels"], title="Atoms")

		fig.canvas.mpl_connect("pick_event", atom_onpick)

		for xs, ys, zs in zip(lines["x"], lines["y"], lines["z"]):
			ax.plot(xs, ys, zs, c="black", lw=1)
	
		ax.set_xlabel("X")
		ax.set_ylabel("Y")
		ax.set_zlabel("Z")
		plt.title("{}".format(self.file_name))
		plt.legend(leg["handlers"], leg["labels"], "upper right")
		if file is not None:
			plt.savefig(file, bbox=0.0)
		else:
			plt.show()

	def load_atoms(self):
		if not self.data:
			if self.verbose:
				print "Loading atoms for empty data"

			return

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
	
	def get_internal_loops(self):
		"""Internal loops defined as any element N, surrounded by C's"""
		loops = []
		if self.verbose:
			print "Looking for internal loops"

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
				self.internal_loops.append(atom)

		if self.verbose:
			print "Found {} internal loops".format(len(loops))
			#pprint(loops)

		self.get_branches()

	def add_branch_from_atom(self, atom, parent):
		if atom in self.internal_loops:
			# create new branches
			# a branch is acomposed of an internal loop id and a neighbour id

			for bond in atom.bonds:
				if bond == parent:
					continue
				if self.verbose:
					print "Found a branch at {}, connected to {}".format(atom.seq_id, bond)
				self.branches.append([atom.seq_id, bond])


	def get_branches(self):
		if not self.internal_loops:
			return
		if self.verbose:
			print "Computing branches"
		# Start at an end node
		# Traverse until internal loop
		
		prev_atom = self.find_first()
		curr_atom = self.find_next(prev_atom)
		while True:
			try:
				parent_id = prev_atom.seq_id
				prev_atom = curr_atom
				curr_atom= self.find_next(curr_atom, parent_id)
			except StopIteration:
				if self.verbose:
					print "Found {} branches".format(len(self.branches))
				break
			else:
				self.add_branch_from_atom(curr_atom, prev_atom.seq_id)
		
		

	def rotate_branch(self, branch, R):
		"""Rotates a branch by a specified rotation matrix"""
		# Get the branch
		branch = self.branches[branch]
		
		loop_atom = self.get_atom_by_id(branch[0])
		loop_child = self.get_atom_by_id(branch[1])
		
		rotation_axis = loop_child.coordinates - loop_atom.coordinates
		
		# Get the subtree
		st = self.get_subtree(loop_atom, loop_child)

		# Rotate all atoms of all sub branches
		for atom in st:
			atom.rotate(R, loop_atom.coordinates, rotation_axis)

		
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


	def find_end_atoms(self):
		return [atom for atom in self.atoms if len(atom.bonds)==1]


	def find_first(self):
		"""Look for node with only 1 bond"""

		for atom in self.atoms:
			if len(atom.bonds) == 1:
				return atom
	
	def find_next(self, atom, parent=None):
		
		for bond in atom.bonds:
			if parent == bond:
				continue
			return self.get_atom_by_id(bond)

		raise StopIteration






	def find_next_loop(self, atom, parent=None):
		"""skip through the molecule bonds until an internal loop is found"""
		for bond_id in atom.bonds:
			if parent and parent == bond_id:
				continue
			bond = self.get_atom_by_id(bond_id)

			if bond in self.internal_loops:
				return bond, atom.seq_id

			return self.find_next_loop(bond, parent=atom.seq_id)

		return None, None

	def get_atom_by_id(self, atom_id):
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
			
		if not os.path.isfile(self.file_path):
			raise ValueError("File {} not found".format(self.file_path))
		self.data = PDB.PDBFile()
		self.data.load_file(open(self.file_path, 'r'))
		self.load_atoms()
		self.get_internal_loops()
		if self.verbose:
			print "Loaded file {}".format(self.file_path)
	

	def write(self, file_name):
		self.update_data()
		if self.verbose:
			print "Writing to file {}".format(file_name)

		self.data.save_file(open(file_name, 'w'))


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
		# Normalize 
		rc_norm = rc/numpy.sqrt(numpy.dot(rc, rc))
		


		# Coordinates must be rotated to be in line with the parent
		# This rotation is along the vector perpendicular to both rc and direction
		rot_unit = numpy.cross(rc_norm, direction)
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
	"""Euler-Rodriguez angles"""
	axis = axis/numpy.sqrt(numpy.dot(axis,axis))

	a = numpy.cos(theta/2)
	b,c,d = -axis*numpy.sin(theta/2)
	return numpy.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
		

def rotation_permutations_from_file(pdb_file, rotation="cubic", verbose=False, plot=False, out_dir=None, interactive=False):


	print "Running rotational permutations on file:  {}".format(pdb_file)
	
	# Read file and find internal loops

	pdb_orig = PDBfile(pdb_file, verbose=verbose)

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

	# Get the rotations for this round
	if rotation == "triangular":
		Rs = rotation_matrices(N=8)
	else:
		Rs = rotation_matrices(N=4)

	from itertools import permutations, combinations_with_replacement

	# For all unique branches from internal loops
	# Compute permutations with the R rotation matrices at loop
	branches = pdb_orig.branches
	if verbose:
		print "A total of {} branches to permute".format(len(branches))
	
	# returns all branch rotation permutations
	#rot_perms = permutations(xrange(len(Rs)), len(branches))
	rot_perms = combinations_with_replacement(xrange(len(Rs)), len(branches))
	for rot_perm in rot_perms:
		f = copy.deepcopy(pdb_orig)
		rot_name = "_".join(["I{}R{}".format(N, R) for N, R  in enumerate(rot_perm)])
		file_name = base_name + rot_name
		if verbose:
			print "Rotation {}".format(rot_name)
		
		for N, R in enumerate(rot_perm):
			# Apply rotation R to internal loop branch N
			f.rotate_branch(N, Rs[R])

		# write the file and plot it
		f.write(os.path.join(file_out_dir, file_name+".pdb"))
		if plot:
			if interactive:
				plot_file = None
			else:
				plot_file = os.path.join(plot_out_dir, file_name + plot_extension)
			f.plot_molecule(file=plot_file)
		del f

	if verbose:
		print "Completed permutations on file {}".format(pdb_file)	


def tests():
	basedir=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
	test_file = os.path.join(basedir, "data", "initialGraphs", "1B36_A_Graph.pdb")
	#test_file = os.path.join(basedir, "data", "initialGraphs", "1GID_A_Graph.pdb")
	print "Using file {}".format(test_file)	
	#test_pdb = PDBfile(test_file, verbose=True)
	#pprint(test_pdb.atoms)
	
	#print test_pdb
	#test_pdb.internal_loops
	#rots = rotation_matrices(N=4, verbose=True)
	#pprint(len(rots))
	#pprint(rots)
	#rotation = "cubic"

	#test_pdb.write(os.path.join(basedir, "tmp", "test.pdb"))
	
	#first_loop = test_pdb.internal_loops[0]
	
	#parent = test_pdb.get_atom_by_id(first_loop.bonds[1])
	#print first_loop
	#print parent
	#subtree = test_pdb.get_subtree(first_loop, parent)
	#pprint(subtree)

	Rs = rotation_matrices(N=4, verbose=True)
	pprint(Rs)
	#R = Rs[0]

	#tt = copy.deepcopy(test_pdb)
	#first_a = tt.atoms[0]

	#print first_a
	#first_a.rotate(R)
	#print first_a

	#test_pdb.plot_molecule(file='./tmp/test.molecule.pdf')
	#test_pdb.plot_molecule()

	out_dir = "tmp/out"
	if os.path.isdir(out_dir):
		os.system("rm -R {}".format(out_dir))
		os.makedirs(out_dir)


	rotation_permutations_from_file(test_file, rotation="triangular", verbose=True, out_dir=out_dir,  plot=True, interactive=True)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Program to parse a PDB file, identify isolation loops, and permute molecular rotations around those loops and write back to a set of PDB files")
	parser.add_argument("-v", "--verbose", help="Print details to console", action="store_true")
	parser.add_argument("-d", "--directory", help="Parse all PDB files in this directory", action="store")
	parser.add_argument("-o", "--output", help="Output new PDB files to this directory", default=os.getcwd(), action="store")
	parser.add_argument("-f", "--file", help="Single PDB file to process", action="store")
	parser.add_argument("-c", "--cubic", help="Rotate around a cubic structure, i.e. 90 deg", action="store_true")
	parser.add_argument("-t", "--triangular", help="Roatate around a triangular structure, i.e. 45 deg", action="store_true")
	parser.add_argument("-p", "--plot", help="Plot the rotated molecules in a `plots` subfolder", default=False, action="store_true")
	parser.add_argument("-i", "--interactive", help="Plot figures interactively", action="store_true")
	parser.add_argument("--test", help="Run the testing script and exit", action="store_true", default=False)
	args = parser.parse_args()

	print
	print "Permuting Rotations v{}".format(__version__)
	print 
	
	if args.test:
		
		print "Running tests"
		tests()
		sys.exit(0)
		
	if args.verbose:
		print "Verbose enabled"

		print "Arguments received: {}".format( args)

	print "Saving output files to directory {}".format(args.output)	
	if not os.path.isdir(args.output):
		if args.verbose:
			print "Creating {}".format(args.output)
		os.mkdir(args.output)

	rotation_method = "cubic"
	if args.triangular:
		rotation_method = "triangular"
	if args.verbose:
		print "Using {} rotation structure".format(rotation_method)

	if args.interactive:
		args.plot = True

	if args.directory:
		if not os.path.isdir(args.directory):
			print "Uh oh, the directory {} does not exist. Exiting".format(args.directory)
			sys.exit(1)

		print "Scanning directory {} for .pdb files".format(args.directory)

		for f in os.listdir(args.directory):
			pdb_file = 	os.path.join(args.directory, f)
			rotation_permutations_from_file(pdb_file, rotation=rotation_method, verbose=args.verbose, plot=args.plot, out_dir=args.output, interactive=args.interactive)


			# Run perms on file

	elif args.file:
		if not os.path.isfile(args.file):
			print "File {} does not exist. Exiting"
			sys.exit(1)
		
		rotation_permutations_from_file(f, rotation=rotation_method, verbose=args.verbose, plot=args.plot, out_dir=args.output, interactive=args.interactive)
		#print "Running rotational permutations on file: {}".format(args.file)
			
	print "All.  Have a nice day"


