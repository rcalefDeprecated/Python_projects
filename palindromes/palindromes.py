#! /usr/bin/env python
"""This module contains several functions useful for analysis of DNA palindromes.

palindromes.py
BME 205 Fall 2014, Programming Assignment #4
November 7th, 2014
Robert Calef

This module contains four functions intended for use in analysis of DNA palindromes.
For many models of palindromes, kmer count data is often needed, typically for a 
range of values for k. The functions provided in this module facilitate the gathering
of this data, and are listed below:

reverse_complement:      Simply takes in a DNA sequence as a string, and returns the 
                         reverse complement of the input sequence.
get_all_dna_palindromes: This function yields all possible DNA palindromes of a given
                         length range, one palindrome at a time.
kmers_from_sequence:     Given a sequence and a value of k, yields each kmer from the
                         sequence, one at a time.
count_kmer_range:        Given a sequence and a range of values for k, returns a Counter
                         object containing counts of occurences of all kmers in the 
                         specified length range in the input sequence.

For more detailed information on each of the functions, see their individual docstrings.
"""
from __future__ import print_function
import sys
import string
from itertools import product
from collections import Counter

#Define a complement table for the following function to speed up reverse complementing
complement_table= string.maketrans("ACGT","TGCA")
def reverse_complement(dna_seq):
    """Takes in a DNA sequence with canonical bases (ACGT) and returns its reverse complement.

    Input:
        dna_seq - A string containing the DNA sequence to be reverse complemented
                  only canonical bases (ACGT) will be complemented, degenerate
                  nucleotide symbols will remain in the reversed sequence as is.
    Output:
        rev_comp - A string containing the reverse complement of 'dna_seq'.

    reverse_complement is a simple function that takes in a DNA sequence as a string,
    and returns the reverse complement of the sequence.
    """
    return dna_seq[::-1].translate(complement_table)

def get_all_dna_palindromes(max_length,min_length=2):
    """A generator that yields all DNA palindromes in a given length range, including odd palindromes.

    Input:
        max_length - An integer specifying the maximum length palindromes to 
                     be yielded.
        min_length - An integer specifying the minimum length palindromes to
                     be yielded, if not specified defaults to 2.
    Output:
        palindrome - A DNA palindrome with length within the specified range.
                     An even length palindrome is defined as some sequence
                                        ABCC'B'A' 
                     where X' indicates the complement of X. We also allow
                     for odd palindromes which are of the form:
                                        ABCXC'B'A'
                     where X is one of the two wildcards specifying base
                     pairing nucleotides, W for A or T, and S for C or G.

    get_all_dna_palindromes is a generator that yields all possible DNA palindromes
    in a given length range. We allow for odd length palindromes, that are simply
    length-1 palindromes with one of the four canonical bases ACGT inserted in the
    center of the palindrome. I.E. we allow the following as valid palindromes:
      FASTA sequence:      5' AATATT   ...   AATCATT   ...  AATAATT 3'
      Reverse complement:  3' TTATAA   ...   AATGATT   ...  TTATTAA 5'
    
      yielded as:             AATATT         AATSATT        AATWATT

    """
    #Negative lengths make no sense in this context.
    if min_length < 0 or max_length < 0:
        print("ERROR: Length cannot be negative.",file=sys.stderr)
        sys.exit(1)
    #Each palindrome consists of length/2 bases, followed by its reverse complement
    #if even, or followed by one of two possible base pairs, and then the reverse
    #complement if odd, hence we generate each possible combination of "ACGT" for
    #length/2, and use this to construct palindromes.
    for length in xrange(min_length,max_length+1):
        for kmer in product("ACGT",repeat=length/2):
            #Need to use join as product returns a tuple, not a string.
            kmer=''.join(kmer)
            if length % 2 == 1:
                for base in "WS":
                    palindrome = kmer + base + reverse_complement(kmer)
                    yield palindrome
            else:
                palindrome = kmer + reverse_complement(kmer)
                yield palindrome


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
    #last kmer is at position n - (k+1) where n is the length of the sequence.
    for start in xrange(len(sequence) - k + 1):
    #Each kmer runs from its start position to k + the start position.
            yield sequence[start:start+k]

def count_kmer_range(sequence,k_min,k_max,ignore_case=True):
    """Takes in a sequence and start and stop characters, and returns counts of kmers in a Counter object.

    Inputs:
        sequence    - A string containing the sequence to be split into kmers,
                      without any start or stop characters.
        k_min       - The minimum length of kmers to be counted, i.e. k_min = 2
                      means that individual character frequencies will not be 
                      counted.
        k_max       - The maximum length of kmers to be counted, i.e. k_max = 4
                      means that 4-mers are the largest kmers that will be 
                      counted.
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
    #counts will store each kmer's count of occurences i.e. counts["ABC"] will be
    #the count of occurences of the 3-mer "ABC" in 'sequence'.
    counts = Counter()
    for k in xrange(k_min,k_max+1):
        for kmer in kmers_from_sequence(sequence,k):
            if ignore_case: kmer = string.upper(kmer)
            counts[kmer] += 1 
    return counts

