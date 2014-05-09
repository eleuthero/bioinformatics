#!/usr/bin/python

from Bio import Entrez
from Bio import SeqIO

# It is requested that when you use the Entrez database, you provide
# your email address so that they acquire usage statistics and can
# identify users causing heavy traffic.

Entrez.email = "burkhac@students.wwu.edu"

# Use the Nucleotide database (there are at least 30 others accessible
# via Entrez) and specify a search query containing the organism name
# and the gene name.

search_handle = Entrez.esearch(
    db="nucleotide",                        # Use nucleotide (nt) database.
    term="Opuntia[Orgn] AND rpl16[Gene]",   # Here's our boolean search.
    usehistory="y"                          # Behave nicely.
    )

search_results = Entrez.read(search_handle)
search_handle.close()

if search_results:

    if search_results['Count'] > 0:

        print "I found {0} entries in Entrez:".format( search_results['Count'] )
        print

        # Save the list of GI numbers.  We need to convert those to
        # accession numbers.

        idlist = search_results['IdList']

        # e-fetch the search results.

        print "Fetching results...",

        fetch_handle = Entrez.efetch(
            db="nuccore",                         # Use the nuccore database.
            id=idlist,                            # Here's the list of GI numbers.
            rettype="gb",                         # Get GenBank (gb) information.
            retmode="text",                       # Return results as text.
            webenv=search_results["WebEnv"],      # Behave nicely.
            query_key=search_results["QueryKey"]  # Behave nicely.
            )

        print "done."

        if fetch_handle:

            # Organize the result into records, parsing as GenBank (gb):

            records = SeqIO.parse(fetch_handle, "gb")

            for record in records:
                print "{0} : {1}".format( record.id, record.description )

            fetch_handle.close()
    else:

        print "I'm sorry, Entrez returned no results for your query."
