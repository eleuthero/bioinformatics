These are some deliberately simple implementations of Markov chains.  We train a Markov model with a source text, in this case the first five chapters of Jane Austen's _Pride and Prejudice_.  From that training, we generate either a collection of fake words, or an excerpt of fake text.

The quality of the text and words so generated is generally poor, since the Markov model we've implemented is the simplest, or "first-order", model, where the next letter/word we choose for our fake word/text generation depends only on:

- the most recent letter/word we generated, and
- the probabilities of transitions from that letter/word in the source text.

Nevertheless, despite how simple the model is, the fake words tend to be strangely pronounceable and the fake text occasionally generates passages that have a definite, albeit disordered, sensibility.
