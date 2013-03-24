#!/usr/bin/env python
"""
PDB support files gracefully borrored from mmLib <http://pymmlib.sourceforge.net/>
"""
import sys
import os
import copy
import math
from pprint import pprint

try:
	import numpy
except ImportError:
	print "This program requires numpy to be installed.  Try `pip install numpy`"
	sys.exit(1)

import PDB
import rotations

__version__ = 0.1
precision = 2 	# Number of decimals to save in float values
plot_extension = ".pdf"

def Rz(alpha):
	"""Rotation matrix around the z axis"""

	cosa = math.cos(alpha)
	sina = math.sin(alpha)

	return numpy.array([[cosa, -sina, 0], [sina, cosa, 0],[0, 0, 1]])

def Ry(alpha):
	"""Rotation matrix around the y axis"""
	cosa = math.cos(alpha)
	sina = math.sin(alpha)

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

	
	matrices = []

	
	for ny in xrange(-N/4, N/4+1):
		if (ny == -N/4) or (ny == N/4):
			if verbose:
				print "Ny: {}".format(ny)
			matrices.append(Ry(ny*step))
			continue
		for nz in xrange(N):
			if verbose:
				print "Nz: {} / Ny: {}".format(nz,ny)

			if (nz*step == numpy.pi/2) and (ny == 0 ):
				# Ignore the rotation that turns it on the incoming bond
				if verbose:
					print "Ignoring nz {} / ny {}".format(nz, ny)
				continue

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
		
		# plot carbons
		carbon_atoms = [atom for atom in self.atoms if atom.element=="C"]
		nitrogen_atoms = [atom for atom in self.atoms if atom.element == "N"]
		other_atoms = [atom for atom in self.atoms if (atom not in carbon_atoms) or (atom not in nitrogen_atoms)]

		# get axes
		fig = plt.figure()
		ax = fig.add_subplot(111, projection="3d")

		carbon = {"colour": "black",
				"size": 400}
		nitrogen = {"colour": "red",
				"size": 320}
		other = {"colour": "grey",
				"size":180}

		in_loop = {"colour": "w",
				"size": 180}
		x = []
		y = []
		z = []
		c = []
		s = []
		lines = {"x":[], "y":[], "z":[] }

		for atom in self.atoms:
			x.append(atom.X)
			y.append(atom.Y)
			z.append(atom.Z)

			if atom in self.internal_loops:
				c.append(in_loop["colour"])
				s.append(in_loop["size"])
			elif atom.element == "C":
				c.append(carbon["colour"])
				s.append(carbon["size"])
			elif atom.element == "N":
				c.append(nitrogen["colour"])
				s.append(nitrogen["size"])
			else:
				c.append(other["colour"])
				s.append(other["size"])
			
			for bond_id in atom.bonds:
				b_atom = self.get_atom_by_id(bond_id)
				lines["x"].append([atom.X, b_atom.X])
				lines["y"].append([atom.Y, b_atom.Y])
				lines["z"].append([atom.Z, b_atom.Z])

		
		#ax.plot(lines["x"], lines["y"], zs=lines["z"], lw=2, c="black")
		#ax.plot(x,y,z, lw=2, c="black")
		ax.scatter(x, y, zs=z, c=c, s=s, marker="o")
		
		for xs, ys, zs in zip(lines["x"], lines["y"], lines["z"]):
			ax.plot(xs, ys, zs, c="black", lw=2)

		plt.title("{}".format(self.file_name))
		#plt.legend(("Carbon", "Nitrogen", "Other"), "upper right")
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
				loops.append(atom)

		
		if self.verbose:
			print "Found {} internal loops".format(len(loops))
			pprint(loops)
		return loops

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

	def find_first(self):
		"""Look for node with only 1 bond"""

		for atom in self.atoms:
			if len(atom.bonds) == 1:
				return atom
	
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

	def get_subtree(self, atom, parent, data=[]):
		"""Returns the molecule after a given atom, traversing away from a parent """
		
		if len(atom.bonds) == 1:
			# End of a chain, return the data tree obtained
			return data
		for bond_id in atom.bonds:
			bonded_atom = self.get_atom_by_id(bond_id)
			if bonded_atom.seq_id == parent:
				continue

			data.append(bonded_atom)
			self.get_subtree(bonded_atom, atom.seq_id, data=data)

		return data



	def add_bond(self, bond):
		"""Add a bond connection"""
		n1 = bond.get("serial")
		n2 = bond.get("serialBond1")
		[atom.bonds.append(n2) for atom in self.atoms if atom.seq_id == n1]
		[atom.bonds.append(n1) for atom in self.atoms if atom.seq_id == n2]

	def read(self):
			
		if not os.path.isfile(self.file_path):
			raise FileError("File {} not found".format(self.file_path))
		self.data = PDB.PDBFile()
		self.data.load_file(open(self.file_path, 'r'))
		self.load_atoms()
		self.internal_loops = self.get_internal_loops()
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

		rc = self.coordinates - point

		# Shift coordinates having the direction vector point along the x-axis
		R  = numpy.identity(3)	
		rotations.R_2vect(R, numpy.array([1,0,0]), direction)
		rc = numpy.dot(R, rc)
		
		# apply rotation matrix to this rotated point
		rc = numpy.dot(matrix, rc)

		# Revert back to original coordinates
		rc = numpy.dot(numpy.linalg.inv(R), rc)
		rc = rc + point

		# Save coordinates to atom
		self.X = round(rc[0], precision)
		self.Y = round(rc[1], precision)
		self.Z = round(rc[2], precision)

		

		

def rotation_permutations_from_file(pdb_file, rotation="cubic", verbose=False, plot=False, out_dir=None):


	print "Running rotational permutations on file:  {}".format(pdb_file)
	
	# Read file and find internal loops

	f = PDBfile(pdb_file, verbose=verbose)

	# If there are loops, then prepare the sub directory to save 
	if len(f.internal_loops) > 0:
		base_name = os.path.splitext(f.file_name)[0]

	# Check the output directory
	if not out_dir:
		out_dir = os.path.dirname(pdb_file)
	file_out_dir = os.path.join(out_dir, base_name)

	if not os.path.isdir(file_out_dir):
		if verbose:
			print "Output Directory not found.  Creating: {}".format(file_out_dir)
		os.mkdir(file_out_dir)
	
	plot_out_dir = os.path.join(file_out_dir, "plots")
	if not os.path.isdir(plot_out_dir) and plots:
		os.mkdir(plot_out_dir)

	# Get the rotations for this round
	if rotation == "triangular":
		Rs = rotation_matrices(N=8)
	else:
		Rs = rotation_matrices(N=4)

	# Start at an end atom
	# Find first internal_loop from this end atom.  
	# Calculate the parent bond direction.  
	# For each rotation matrix, 
	# 	copy the molecule, 
	# 	rename it accordingly
	#	apply the rotation at the internal loop to the subtree
	# 	recurse into the next internal loop for each rotation
	#	if an end atom is found:
	#		write the current molecule to PDB file
	#		(Optional) Plot the 3d structure
	
	print "Applying rotation permutations"
	iterate_rotations(f, Rs, output_dir=file_out_dir, plot_dir=plot_out_dir, plot=plot)

def iterate_rotations(pdb_file, Rs, il_atom=None, il_parent=None, il_count=0, output_dir=None, plot_dir=None, plot=False, verbose=False):
	base_name = os.path.splitext(pdb_file.file_name)[0]
	if not il_atom:
		il_atom, il_parent = pdb_file.find_next_loop(pdb_file.find_first())
		il_count += 1
		if verbose:
			print
			print "Starting iterations at loop atom: {} from parent ".format(il_atom, il_parent)

	for r, R in enumerate(Rs):

		f_name ="{}_I{}_R{}".format(base_name, il_count, r)
		
		if verbose:
			print "Rotation Iteration: {}".format(f_name)

		# copy the tree
		fc = copy.deepcopy(pdb_file)
		
		# Apply the rotation to the copied tree
		subtree = fc.get_subtree(il_atom, il_parent)
		for atom in subtree:
			fc_atom = fc.get_atom_by_id(atom.seq_id)
			direction = fc_atom.coordinates - il_atom.coordinates
			fc_atom.rotate(R, il_atom.coordinates, direction)

		# find next il
		next_il_atom, next_il_parent = fc.find_next_loop(il_atom, parent=il_parent)
		if not next_il_atom:

			# we're done, write the tree to a file and/or plot it
			fc.write(os.path.join(output_dir, f_name+".pdb"))
			if plot:
				fc.plot_molecule(file=os.path.join(plot_dir, f_name+plot_extension))

		else:
			#il_count += 1
			iterate_rotations(fc, Rs, il_atom=next_il_atom, il_parent=next_il_parent, il_count=il_counti+1, out_dir=out_dir, plot_dir=plot_dir, plot=plot)

def tests():
	basedir=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
	test_file = os.path.join(basedir, "data", "initialGraphs", "1B36_A_Graph.pdb")

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

	#Rs = rotation_matrices(N=4)
	#R = Rs[0]

	#tt = copy.deepcopy(test_pdb)
	#first_a = tt.atoms[0]

	#print first_a
	#first_a.rotate(R)
	#print first_a

	#test_pdb.plot_molecule(file='./tmp/test.molecule.pdf')

	rotation_permutations_from_file(test_file, verbose=True, out_dir="tmp/out", plot=True)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Program to parse a PDB file, identify isolation loops, and permute molecular rotations around those loops and write back to a set of PDB files")
	parser.add_argument("-v", "--verbose", help="Print details to console", action="store_true")
	parser.add_argument("-d", "--directory", help="Parse all PDB files in this directory", action="store")
	parser.add_argument("-o", "--output", help="Output new PDB files to this directory", action="store")
	parser.add_argument("-f", "--file", help="Single PDB file to process", action="store")
	parser.add_argument("-c", "--cubic", help="Rotate around a cubic structure, i.e. 90 deg", action="store_true")
	parser.add_argument("-t", "--triangular", help="Roatate around a triangular structure, i.e. 45 deg", action="store_true")
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

	if args.output:
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

	if args.directory:
		if not os.path.isdir(args.directory):
			print "Uh oh, the directory {} does not exist. Exiting".format(args.directory)
			sys.exit(1)

		print "Scanning directory {} for .pdb files".format(args.directory)

		for f in os.listdir(args.directory):
			
			rotation_permutationss_from_file(f, rotation=rotation_method, verbose=args.verbose)


			# Run perms on file

	if args.file:
		if not os.path.isfile(args.file):
			print "File {} does not exist. Exiting"
			sys.exit(1)
		
		rotation_permutations_from_file(f, rotation=rotation_method, verbose=args.verbose)
		#print "Running rotational permutations on file: {}".format(args.file)
			
	print "Finished.  Have a nice day"


