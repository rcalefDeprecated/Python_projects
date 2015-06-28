#! /usr/bin/env python2.7
"""This module contains functions useful for constructing a Markov model for FASTA sequences.

markov.py
BME 205 Fall 2014, Programming Assignment #3
October 31st, 2014
Robert Calef

This module contains four functions intended for use in constructing an order k 
stochastic Markov model for sequences of characters. In order to estimate parameters 
for an order k Markov model in a simple manner, one typically wants to count the
occurences of k-mers over some training set of sequences. The functions provided in
this module can be used to gather this data, and are listed below:

count_kmers:         Takes in a sequence of characters, and returns a Counter object 
                 containing counts of occurences of all kmers in the sequence.
kmers_from_sequence: Splits a sequence into kmers and yields each kmer in the sequence, one 
                     at a time.
generate_all_kmers:  Takes in an alphabet, a start character, and a stop character,
                     and yields all valid kmers, one at a time.
check_alphabet:      Checks an alphabet, and start and stop characters to make sure they 
                     are valid.

For more information on each of the functions, see their individual docstrings.
"""
from __future__ import print_function
import string
from collections import Counter
from itertools import product

def count_kmers(sequence,k,start="^",stop="$",ignore_case=True):
    """Takes in a sequence, and start and stop characters, and returns counts of kmers in a Counter object.

    Inputs:
        sequence    - A string containing the sequence to be split into kmers,
                      without any start or stop characters.
        k           - The length of kmers to be counted, i.e. k = 1 will count 
                      occurences of individual characters,
        start       - The character used to represent the start of the sequence.
                      Defaults to '^'.
        stop        - A character representing the end of the sequence.
                      Defaults to '$'.
        ignore_case - A Boolean specifying whether or not to preserve case of 
                      letters in 'sequence', defaults to True. If set to False,
                      the 3-mers "AAA" and "AaA" will be treated as two distinct
                      3-mers, if True, all kmers will be converted to uppercase.
    Output:
        counts   - A Counter object containing (kmer, counts of occurences) pairs
                   as (key, value) pairs. That is, counts["AAA"] will return the
                   number of times the kmer "AAA" occurs in 'sequence'.
    
    count_kmers will split 'sequence' in to kmers of length 'k', prepending and 
    appending the start and stop characters as necessary, and the occurences of
    each kmer will be counted. Only kmers that occur in the sequence will have 
    an entry in the returned Counter object.
    """
    #counts["AAA"] will be the number of times the 3-mer "AAA" was seen in 'sequence'
    counts = Counter()
    if k < 3: 
        #If k=1 or k=2, then we only want one start character and one stop character,
        #as we don't want kmers consisting of only start characters or only stop characters
        #unless k=1, in which case we count each character, including start and stop
        #characters individually.
        if ignore_case: seq = string.upper(start + sequence + stop)
        else: seq = start + sequence + stop
    else:
        #If k >2, then we want to prepend and append k-1 start and stop characters to
        #the sequence, so the first kmer will consist of k-1 start characters followed
        #by the first character of the sequence, and the last kmer will be the final
        #character of the sequence followed by k-1 stop characters.
        if ignore_case: seq = string.upper(start*(k-1) + sequence + stop*(k-1))
        else: seq = start*(k-1) + sequence + stop*(k-1)
    #Finally, we just get each kmer from the sequence, and increment the kmer's count
    for kmer in kmers_from_sequence(seq,k):
        counts[kmer] += 1
    return counts


def kmers_from_sequence(sequence,k):
    """A simple generator to split a sequence in to kmers, yielding one kmer at a time.
    
    Inputs:
        sequence - A string containing the sequence to be split into kmers.
        k        - An integer specifying the length of kmers.
    Outputs:
        kmer     - A length 'k' substring of 'sequence'.

    kmers_from_sequence is a short generator that takes a string 'sequence', 
    and an integer 'k' as input, and yields each length 'k' substring of 'sequence'
    one at a time.
    """
    #The start of the first kmer is position 0 in the sequence, and the start of the 
    #last kmer is at position n - (k+1) where n is the length of the sequence. Each
    #kmer runs from its start position to k + the start position.
    for start in xrange(len(sequence) - k + 1):
        yield sequence[start:start+k]

def generate_all_kmers(k,alphabet,start,stop):
    """Generates all possible kmers over the given alphabet, including start and stop characters.
    
    Inputs:
        k        - An integer specifying the length of kmers to be generated.
        alphabet - Any iterable containing the alphabet of valid characters
                   for kmers, should not contain the start and stop characters.
        start    - The character used to represent the start of the sequence.
        stop     - The character used to represent the end of the sequence,
    Outputs:
        kmer     - A string of length 'k' consisting of characters from the
                   specified alphabet, and the start and stop characters.
               All such valid kmers will be generated.
    
    generate_all_kmers takes in an alphabet, start and stop characters, and an 
    integer 'k', and yields all valid kmers over this alphabet, one kmer at a 
    time. A valid kmer comes from one of the following four cases (in the descriptions
    below, n is an integer such that 0<n<k). Examples are given for an alphabet of
    "ABC", '^' as the start character, '$' as the stop character, and k=4:
        1. A string of length k consisting only of characters from the alphabet.
            e.g. "AAAA", "AAAB", "AAAC" ...
        2. k-n characters from the alphabet prefixed by n start characters.
            e.g. "^AAB", "^^AA", "^^^A", BUT NOT "^^^^"
        3. k-n characters from the alphabet suffixed by n stop characters.
            e.g. "CCA$", "CA$$", "A$$$", BUT NOT "$$$$" 
        4. k-n characters from the alphabet prefixed by x start characters
           and suffixed by y stop characters, where x + y = n, and x and y are
           both at least 1. For this case only, n is allowed to equal k, so as 
           to generate kmers corresponding to an empty sequence.
                e.g. "^AA$", "^AB$", "^^^$", "^^$$", "^$$$"
    """
    #First we yield all possible strings of length k with characters drawn from
    #the specified alphabet. The itertools.product function is used to generate
    #all of our kmers. itertools.product(alphabet,repeat=k) takes in our alphabet,
    #and yields all elements of the Cartesian product of the alphabet with itself 
    #'k' times. product() yields tuples, so we use ''.join() to get kmers as strings.
    for kmer in product(alphabet,repeat=k):
        yield ''.join(kmer)
    #Next we want to yield kmers containing k-n characters drawn from the alphabet.
    for num in xrange(1,k):

        #Generate all possible (k-n)-mers, and use these to get all kmers from categories
        #2 and 3 as described in the fuctio docstring.
        for sub_mer in product(alphabet,repeat=k-num):
            sub_mer = ''.join(sub_mer)
            kmer = (num * start) + sub_mer
            yield kmer
            kmer= sub_mer + (num * stop)
            yield kmer

        #For kmers containing both start and stop characters, we will use 'num'
        #as the number of start characters, and 'num_stops', which ranges
        #from 1 to (k - num + 1), as we want at most (k-num) stop characters
        #hence (k-num+1) as the upper bound, since xrange uses a half-open
        #interval. Note that the number of alphabet characters will be 
        #k -  num - nnum_stops
        for num_stops in xrange(1,k-num + 1):
            mid_length = k - num - num_stops
            #If mid_length=0, product returns an empty tuple, so we will
            #yield kmers corresponding to the empty sequence.
            for middle in product(alphabet,repeat=mid_length):
                kmer = (start * num) + ''.join(middle) + (stop*num_stops)
                yield kmer

def check_alphabet(alphabet,start,stop):
    """Checks alphabet for invalid characters, prints error if such characters are found.

    Inputs:
            alphabet - Any iterable supporting the 'in' operator and containing the
                       alphabet of allowed characters for FASTA sequences.
            start    - The character used to represent the start of a FASTA sequence.
            stop     - The character used to represent the end of a FASTA sequence.
    Output:
            none, only prints errors if problems with the alphabet are found.

    check_alphabet() is a small utility function that performs checks to make sure the
    alphabet is valid, that is, that the alphabet contains only printable,
    non-whitespace characters, and does not contain the specified start or stop
    characters. We must make sure the start and stop characters are not in the
    alphabet so as to make sure the start and stop characters do not appear in the
    middle of any of the FASTA sequences. Start and stop characters in the middle
    of sequences would lead to erroneous start and stop k-mer counts.
    """
    if start in alphabet or stop in alphabet:
        print("ERROR: Alphabet contains start and/or stop characters.\n",
            file=sys.stderr)
        sys.exit(1)

    #To check if alphabet is valid, we first make a set of valid characters, and then
    #make sure the alphabet is contained within this set. The '-' operator below
    #makes 'allowed_chars' to be a set of all printable, non_whitespace characters.
    #The '<=' operator on sets checks if the left set is a complete subset of the right
    #set, i.e. that every element of the left set is also an element of the right set
    allowed_chars = set(string.printable)-set(string.whitespace)
    if not set(alphabet) <= allowed_chars:
        print("WARNING: Alphabet contains either whitespace "
              "characters, or non-printable characters.",file=sys.stderr)



