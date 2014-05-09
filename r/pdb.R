#!/usr/bin/R

# Load the Bio3D library that makes it easy for us to parse
# the PDB file's HELIX and SHEET records.

library("bio3d")

# Infile.  A list of protein IDs separated by newlines or spaces.

infile  <- "proteins.txt"

# Outfile.  Will be recreated on every run.

outfile <- "proteins_out.txt"

# Delete and recreate the output file.

if (file.exists(file = outfile))
{
    file.remove(file = outfile)
}

file.create(file = outfile)

# Get a list of protein IDS from infile.

pids <- strsplit( readLines(infile), "[[:space:]]+" )

# Write the header to the outfile.

cat("Protein\t% helix\t% sheet\n", file = outfile, append = TRUE)

# Iterate over every protein ID in infile.

for (pid in pids)
{
    # Get the PDB entry for this protein online.

    pdb <- read.pdb(pid)

    # Print the % alpha helix and % beta sheet
    # values based on the number and lengths of the
    # structures compared to the number of residues
    # in the sequence.

    helixpct <- sum(pdb$helix$end - pdb$helix$start) / length(pdb$seqres) * 100
    sheetpct <- sum(pdb$sheet$end - pdb$sheet$start) / length(pdb$seqres) * 100

    formatted <- sprintf("%s\t%.2f\t%.2f\n", pid, helixpct, sheetpct)

    cat(formatted, file = outfile, sep = "", append = TRUE)
}
