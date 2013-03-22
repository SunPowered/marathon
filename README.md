Rotational Permutations of RNA Folds
=====================================

This is the README

# Goal

The intent of this project is to create a script that is able to parse protein database (_.pdb_) files
and compute the various rotational permutations at determined initial loop sites.  The rotation
can either be in increments of 45 or 90 degrees in all 3 dimensions.  The various folded structures
are then written back to a _.pdb_ format.

# Usage

The program is run from the command line using python 2.7.  It takes the following command line arguments:

	* directory: Parse all .pdb files in the given directory
	* output: Specifies the output directory to write new files
	* file: A single .pdb file to parse
	* verbose: Print more information to the console

## Examples

`rotate_pdb.py --output myoutputdir --file mydatafile.pdb`
	
will output all rotational permuatations of `mydatafile.pdb` to `myoutputdir`



# Help

Run `rotate_pdb.py --help`

