#!/usr/bin/python

import re
import sys

freq = { }

def getCorpus(filename):
    clean_regex = re.compile(r"[^A-Za-z\.\?\!]")
    split_regex = re.compile(r"[.?!]")

    with open(filename) as fin:
        corpus = fin.read()

    # Clean up the corpus by removing most punctuation.
    # Leave the sentence-terminating punctuation; the last
    # letter of a sentence should not be interpreted as having
    # any effect on the first letter of the next sentence.

    corpus = re.sub(clean_regex, '', corpus).lower()

    # Return the sentences in the corpus.

    return re.split(split_regex, corpus);

def compileFrequency(current, next):
    global freq

    # Have we seen the current letter before ?

    if current in freq:
        entry = freq[current]
    else:
        entry = { }
        freq[current] = entry

    # Have we seen the transition from the current to the next letter ?

    if next in entry:
        entry[next] += 1
    else:
        entry[next] = 1

def processSentence(sentence):

    # Only process sentences of at least two letters.

    if len(sentence) > 1:

        # Turn sentence into a queue.

        sequence = list(sentence)
        sequence.reverse()

        current = sequence.pop()

        while len(sequence) > 0:
            next = sequence.pop()

            # Add current --> next to our transition tally.

            compileFrequency(current, next)
            current = next 

def showFrequenciesList():
    global freq

    for current in sorted(freq):
        print current, ' -->'

        sum = 0

        for next in freq[current]:
            sum += freq[current][next]

        for next in sorted(freq[current]):
            count = float(freq[current][next])
            print '\t%s: %5i / %5i (%.3f%%)' % (next, count, sum, count * 100 / sum)

        print

def showFrequenciesTable():
    global freq
    keys = list('abcdefghijklmnopqrstuvwxyz')

    # Print header row.

    print 'NXT CUR,',
    for current in keys:
   	print current, ',',
    print

    print '========',
    for current in keys:
        print '=======,',
    print

    # For each letter in the transition set...

    for next in keys:

        # Print header cell

        print next, ',',

        # For each transition into the next letter ...

        for current in keys:
            if next in freq[current]:

                # Calculate probability of current --> next transition

                sum = 0
                for transition in freq[current]:
                    sum += freq[current][transition]
                print '%.5f%%,' % (100 * float(freq[current][next]) / sum), 
            else:
                print '0,',
        print

# ====
# Main
# ====

# Ensure user has given us a corpus filename.

if len(sys.argv) < 2:
    print "Enter corpus file as first argument to script."
    exit(1)

# Calculate frequency of letter transitions in the corpus.

sentences = getCorpus(sys.argv[1])

for sentence in sentences:
    processSentence(sentence)

# Show transition analysis.

# showFrequenciesList()      # Enable this to see a list of transitions
showFrequenciesTable()       # Enable this to see a table of transitions
