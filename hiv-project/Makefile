

PYTHON=/usr/bin/python

all:
	mkdir ./sequences
	$(PYTHON) ./1_process.py
	$(PYTHON) ./2_consensus.py
	$(PYTHON) ./3_import.py

clean:
	rm -rf ./sequences	
