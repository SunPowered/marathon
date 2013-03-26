.PHONY: all test clean

OUT_DIR=tmp/out
SCRIPT=src/marathon.py
TEST_SCRIPT=test/test_marathon.py

remove:
	rm -Rf $(OUT_DIR)

test: 
	python $(TEST_SCRIPT)

