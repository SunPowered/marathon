#!/usr/bin/env python

import argparse

__version__ = "0.1"

if __name__ == "__main__":

	import os
	import sys

	parser = argparse.ArgumentParser(description="Find the min rmsd from a formatted text file")

	parser.add_argument('file', help="The file to process")
	parser.add_argument('-v', '--verbose', help="Display more information to the console", default=False, action="store_true")
	args = parser.parse_args()

	if not args.file:
		print "File must be given"
		sys.exit()

	if not os.path.isfile(args.file):
		print "File {} does not exist.  Exiting".format(args.file)
		sys.exit()

	print "min_rmsd v{}".format(__version__)
	print "Finding the minimum RMSD for file {}".format(args.file)

	min_rmsd = 999999999999999
	min_rmsd_counter = []

	f = open(args.file, 'r')

	for line in f:
		try:
			l_s= line.split()
			cur_iter = l_s[0]
			cur_rmsd = float(l_s[1])
		except Exception, e:
			print
			print "An error occurred while reading the file.  Is it not formatted properly?"
			print "Error msg: {}".format(e)
			sys.exit()

		if cur_rmsd < min_rmsd:
			if args.verbose:
				print "Found new RMSD min {} at line {}".format(cur_rmsd, line_counter)
			min_rmsd = cur_rmsd
			min_rmsd_counter = [cur_iter]
		elif cur_rmsd == min_rmsd:
			min_rmsd_counter.append(cur_iter)


	print
	print "Minimum rmsd of {} at iteration(s) {}".format(min_rmsd, min_rmsd_counter)
	print "All done.  Have a nice day"
		
