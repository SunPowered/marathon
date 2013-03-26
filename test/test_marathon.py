#!/usr/bin/env python

import sys
import os
import shutil

base_dir = os.path.dirname(__file__)
marathon_script = os.path.abspath(os.path.join(base_dir, os.pardir, "src", "marathon.py"))

test_data_directory = os.path.join(base_dir, 'initialGraphs')
single_data_file = os.path.join(test_data_directory, "1B36_A_Graph.pdb")

output_directory = os.path.join(base_dir, "output")

if os.path.isdir(output_directory):
	# remove this directory
	shutil.rmtree(output_directory)

os.system("python {} -f {} -o {} -t -v -i".format(marathon_script, single_data_file, output_directory))
 


