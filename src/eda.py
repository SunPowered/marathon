#!/usr/bin/env python
import os
from mmLib import PDB

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
test_file = os.path.join(base_dir, 'data', 'initialGraphs', '1B36_A_Graph.pdb')

pfile = PDB.PDBFile()
pfile.load_file(open(test_file, 'r'))



