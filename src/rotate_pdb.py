#!/usr/bin/env python
"""
PDB support files gracefully borrored from mmLib <http://pymmlib.sourceforge.net/>
"""
import sys
import os

import math
import numpy

import PDB

__version__ = 0.1

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

	
	def Rz(n):
		cosa = math.cos(n*step)
		sina = math.sin(n*step)

		return numpy.array([[cosa, -sina, 0], [sina, cosa, 0],[0, 0, 1]])

	def Ry(n):
		cosa = math.cos(n*step)
		sina = math.sin(n*step)

		return numpy.array([[cosa, 0, -sina], [0, 1, 0], [sina, 0, cosa]])

	matrices = []

	
	for ny in xrange(-N/4, N/4+1):
		if (ny == -N/4) or (ny == N/4):
			if verbose:
				print "Ny: {}".format(ny)
			matrices.append(Ry(ny))
			continue
		for nz in xrange(N):
			if verbose:
				print "Nz: {} / Ny: {}".format(nz,ny)

			if (nz*step == numpy.pi/2) and (ny == 0 ):
				# Ignore the rotation that turns it on the incoming bond
				if verbose:
					print "Ignoring nz {} / ny {}".format(nz, ny)
				continue

			matrices.append(numpy.dot(Rz(nz), Ry(ny)))
	return matrices
	
class PDBfile(object):
	"""PDB File to read modify and write"""

	def __init__(self, file_name, verbose=False):
		"""Arguments: 
			file_name:  The file to read"""
		self.file_name = file_name
		self.verbose = verbose
		self.data = []
		self.atoms = []

		self.read()

	def __str__(self):
		return str(self.data)
	
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
		return loops


	def add_bond(self, bond):
		"""Add a bond connection"""
		n1 = bond.get("serial")
		n2 = bond.get("serialBond1")
		[atom.bonds.append(n2) for atom in self.atoms if atom.seq_id == n1]
		[atom.bonds.append(n1) for atom in self.atoms if atom.seq_id == n2]

	def read(self):
			
		if not os.path.isfile(self.file_name):
			raise FileError("File {} not found".format(self.file_name))
		self.data = PDB.PDBFile()
		self.data.load_file(self.file_name)
		self.load_atoms()
		if self.verbose:
			print "Loaded file {}".format(self.file_name)
	

	def write(self, file_name):
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

def roation_permuations_from_file(pdb_file, rotation="cubic", verbose=False):


	print "Running rotational permutations on file:  {}".format(pdb_file)
	
	# Read file and find internal loops

	f = PDBfile(pdb_file, verbose=verbose)
	loops = f.get_internal_loops()

	# If there are loops, then prepare the sub directory to save 
		

def tests():
	from pprint import pprint
	basedir=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
	test_file = os.path.join(basedir, "data", "initialGraphs", "1B36_A_Graph.pdb")

	test_pdb = PDBfile(test_file, verbose=True)
	print test_pdb.atoms
	#print test_pdb
	#print test_pdb.get_internal_loops()
	pprint(len(rotation_matrices(N=4, verbose=True)))

	rotation = "cubic"

	#test_pdb.write(os.path.join(basedir, "tmp", "test.pdb"))

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
		sys.exit(1)
		
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


