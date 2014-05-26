PYTHON=/usr/bin/python
OUTFILE=schema/hiv_data.sql
MAPFILE=alignment_map.html
THRESHOLD=1

all:
	mkdir ./sequences
	$(PYTHON) ./process.py
	$(PYTHON) ./consensus.py
	$(PYTHON) ./reduce.py $(THRESHOLD)
	$(PYTHON) ./import.py > $(OUTFILE)
	$(PYTHON) ./alignment.py > $(MAPFILE)

clean:
	rm -f $(OUTFILE) $(MAPFILE)
	rm -rf ./sequences	
