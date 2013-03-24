OUT_DIR=tmp/out
SCRIPT=src/rotate_pdb.py

remove:
	rm -Rf $(OUT_DIR)

test: remove
	$(SCRIPT) -d data/initialGraphs -o $(OUT_DIR) -p -v	


