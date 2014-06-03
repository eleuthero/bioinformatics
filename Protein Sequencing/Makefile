PYTHON=/usr/bin/python
OUTFILE=proteins/
THRESHOLD=1

all:
    mkdir./proteins
    $(PYTHON) ./import.py > $(OUTFILE)

clean:
    rm -f $(OUTFILE) $(MAPFILE)
    rm -f "./alignment_*.html"
    rm -rf ./proteins
