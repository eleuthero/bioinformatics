PYTHON=/usr/bin/python
OUTFILE=schema/hiv_data.sql
MAPFILE=alignment_map.html

all:
	mkdir ./sequences
	$(PYTHON) ./1_process.py
	$(PYTHON) ./2_consensus.py
	$(PYTHON) ./3_import.py > $(OUTFILE)
	$(PYTHON) ./4_alignment.py > $(MAPFILE)

clean:
	rm -f $(OUTFILE) $(MAPFILE)
	rm -rf ./sequences	
