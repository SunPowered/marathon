Rotational Permutations of RNA Folds
=====================================


# Goal

The intent of this project is to create a script that is able to parse protein database (*.pdb*) files
and compute the various rotational permutations at determined initial loop sites.  The rotation
can either be in increments of 45 or 90 degrees in all 3 dimensions.  The various folded structures
are then written back to a *.pdb* format.

# Usage

The program is run from the command line using python 2.7.  It takes the following command line arguments:

	* directory: Parse all .pdb files in the given directory
	* output: Specifies the output directory to write new files
	* file: A single .pdb file to parse
	* cubic: Use a cubic structure to rotate around (90 deg)
	* triangular: Use a triangular structure to rotate around (45 deg)
	* verbose (Boolean): Print more information to the console
	* plot (Boolean): Save the permutation structure plots to file
	* interactive (Boolean): Send the plots to the default display to view interactively
	

## Examples

`rotate_pdb.py --output myoutputdir --file mydatafile.pdb -c -p -v`
	
will output all rotational permuatations of `mydatafile.pdb` to `myoutputdir`

# Requirements

The following packages must be installed in the python path:

	* numpy <http://www.numpy.org/>
	* matplotlib (Required for plotting) <http://matplotlib.org/>

# Help

Run `rotate_pdb.py --help` to get a list of available options

# Author

This script was writting by Tim van Boxtel - 2013
