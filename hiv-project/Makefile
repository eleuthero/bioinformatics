PYTHON=/usr/bin/python
OUTFILE=schema/hiv_data.sql
THRESHOLD=1

all:
	mkdir ./sequences
	$(PYTHON) ./process.py
	$(PYTHON) ./consensus.py
	$(PYTHON) ./reduce.py $(THRESHOLD)
	$(PYTHON) ./import.py > $(OUTFILE)
	$(PYTHON) ./alignment.py

clean:
	rm -f $(OUTFILE) $(MAPFILE)
	rm -f "./alignment_*.html"
	rm -rf ./sequences	
