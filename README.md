Rotational Permutations of RNA Folds
=====================================


# Goal

The intent of this project is to create a script that is able to parse protein database (*.pdb*) files
and compute the various rotational permutations at determined initial loop sites.  The rotation
can either be in increments of 45 or 90 degrees in all 3 dimensions.  The various folded structures
are then written back to a *.pdb* format.

# Usage

The program is run from the command line using python 2.7.  It takes the following command line arguments:

* output: Specifies the output directory to write new files
* cubic: Use a cubic structure to rotate around (90 deg)
* triangular: Use a triangular structure to rotate around (45 deg)
* detailed:  Use detailed rotation names for each branch, otherwise the iterations are autoincremented 
* verbose (Boolean): Print more information to the console
* plot (Boolean): Save the permutation structure plots to file
* interactive (Boolean): Send the plots to the default display to view interactively
	

## Examples

`python marathon.py --output myoutputdir -p -v mydatafile.pdb`
	
will output all cubic rotational permuatations of `mydatafile.pdb` to `myoutputdir` with plots for each iteration.

`python marathon.py -o myoutput_dir -t my_mol.pdb`

will outoput all rotational	

# Requirements

The following packages must be installed in the python path:

	* numpy <http://www.numpy.org/>
	* matplotlib (Required for plotting) <http://matplotlib.org/>

# Help

Run `marathon.py --help` to get a list of available options

# Author

This script was writting by Tim van Boxtel - 2013
