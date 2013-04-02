.PHONY: all test clean

SCRIPT=src/marathon.py
TEST_SCRIPT=test/test_marathon.py

remove:
	rm -Rf $(OUT_DIR)

test: 
	python $(TEST_SCRIPT)
	
clean-docs:
	rm doc/marathon.pdf

docs: 
	pandoc -o doc/marathon.pdf -s --toc doc/marathon.md
