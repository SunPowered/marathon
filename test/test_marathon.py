#!/usr/bin/env python

import sys
import os
import shutil

base_dir = os.path.dirname(__file__)
marathon_script = os.path.abspath(os.path.join(base_dir, os.pardir, "src", "marathon.py"))

#test_data_directory = os.path.join(base_dir, 'initialGraphs')
test_data_directory = os.path.join(base_dir, 'testGraphs')
#filename = "1B36_A_Graph.pdb"
filename = "1DUQ_AB_Graph.pdb"
#filename = "1GID_A_Graph.pdb"
single_data_file = os.path.join(test_data_directory, filename)

#output_directory = os.path.join(base_dir, "debug")
output_directory = os.path.join(base_dir, 'output')
#if os.path.isdir(output_directory):
	# remove this directory
	#shutil.rmtree(output_directory)

data_args = single_data_file
#data_args = test_data_directory
os.system("python {} -o {} -d -r -v --print-skips {}".format(marathon_script, output_directory, data_args))
 


