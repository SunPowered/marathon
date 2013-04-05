Marathon
==========

# Goal

The intent of this project is to create a script that is able to parse protein database (*.pdb*) files
and compute the various rotational permutations at determined initial loop sites.  The rotation
can either be in increments of 45 or 90 degrees in all 3 dimensions.  The various folded structures
are then written back to a *.pdb* format.

# Usage

	usage: marathon.py [-h] [-v] [-o OUTPUT] [-c] [-t] [-p] [-i] [-r] [-d]
    	               [--print-skips]
    	               args [args ...]

	Program to parse a PDB file, identify isolation loops, and permute molecular
	rotations around those loops and write back to a set of PDB files

	positional arguments:
	  args                  One or more filenames or directories

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --verbose         Print details to console
	  -o OUTPUT, --output OUTPUT
	                        Output new PDB files to this directory
	  -c, --cubic           Rotate around a cubic structure (default), i.e. 90 deg
	  -t, --triangular      Roatate around a triangular structure, i.e. 45 deg
	  -p, --plot            Plot the rotated molecules in a `plots` subfolder
	  -i, --interactive     Plot figures interactively
	  -r, --rmsd            Calculate the root means square distance of each
	                        iteration and save all values to a file
	  -d, --detailed        Save the iteration names with detailed information for
	                        each branch and rotation number, otherwise just use
	                        the iteration counter as a name
	  --print-skips         Print the skipped rotation iteration names to a file
	                        in the output directory

The program is run from the command line using python 2.7.  

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

Run `python marathon.py --help` to get a list of available options


