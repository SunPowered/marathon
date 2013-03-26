OUT_DIR=tmp/out
SCRIPT=src/marathon.py

remove:
	rm -Rf $(OUT_DIR)

test: remove
	$(SCRIPT) -d data/initialGraphs -o $(OUT_DIR) -p -v	


