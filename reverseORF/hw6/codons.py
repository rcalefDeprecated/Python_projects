#! /usr/bin/env python2.7
"""This module contains functions useful for constructing null models of open reading frames.

codons.py
BME 205 Fall 2014, Programming Assignment #6
November 21st, 2014
Robert Calef

This module contains three functions and a class intended for use in constructing 
stochastic models of open reading frames. The module has kmer counting functions which can
be used to estimate parameters for an ORF model using 3-mer counts from a given set of 
sequences. The functions provided in this module can be used to gather this data, and are 
listed below:

reverse_complement:  Takes in a DNA or RNA sequence containing canonical bases ACGT or 
                     ACGU and returns the reverse complement of the sequence as a string.
count_kmers:         Takes in a sequence of characters, and returns a Counter object 
                     containing counts of occurences of all kmers in the sequence.
kmers_from_sequence: Splits a sequence into kmers and yields each kmer in the sequence, one 
                     at a time.

The codons.py module also defines a class used for generating random codons according to
a codon usage table:

codon_generator:     A class defining an object that can be used to generate codons 
                     randomly according to the data in a codon usage table.

For more information on each of the functions, and the codon_generator class, see their
respective docstrings.
"""
from __future__ import print_function
import string
import sys
import urllib2
import re
from collections import Counter
from random import random
from bisect import bisect


#Define complement tables for the following function to speed up reverse complementing.
dna_complement_table = string.maketrans("AGCT","TCGA")
rna_complement_table = string.maketrans("AGCU","UCGA")
def reverse_complement(sequence):
    """Takes in a DNA or RNA sequence with canonical bases (ACGT or ACGU) and returns its reverse complement.

    Input:
        sequence - A string containing sequence to be reverse complemented
                   only canonical bases (ACGT or ACGU) will be complemented,
                   degenerate nucleotide symbols will remain in the reversed
                   sequence as is. The sequence is considered an RNA sequence
                   if it contains a 'U' in the sequence.
    Output:
        rev_comp - A string containing the reverse complement of 'sequence'.

    reverse_complement is a simple function that takes in a DNA or RNA sequence as a 
    string, and returns the appropriate reverse complement of the sequence. 
    WARNING: This function will only reverse complement canonical bases, ACGT or ACGU
    for DNA and RNA respectively. All other characters will be left unchanged.
    """
    if 'U' in string.upper(sequence):
        return sequence[::-1].translate(rna_complement_table)
    else:
        return sequence[::-1].translate(dna_complement_table)


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




class codon_generator:
    """An object that reads in a codon usage table, and uses this to generate random codons with the appropriate weights.

    The codon_generator class defines an object used for generating random codons 
    according to the codon biases of a specific organism. Codon biases are determined
    from a codon usage table in the format provided by the Codon Usage Database 
    (www.kazusa.or.jp/codon/). In brief, the format consists of one entry per codon, 
    with each entry in the following form:
            codon amino_acid fraction frequency_per_thousand (number)
       
        codon - The codon to which the other data pertains
        amino_acid - The amino acid specified by 'codon' in this organism
        fraction - The fraction of codons for 'amino_acid' made up by 'codon'
        per_thousand - Expected number of occurences of 'codon' in 1000 random codons,
                       typically not an integer.
        number       - The raw count of occurences of 'codon' across all coding sequences
                       used to obtain the codon counts.
    For example, a single entry would look like this:
                                   UUU F 0.59 26.3 ( 22246)      
    While the tables provided by the Codon Usage Database have four entries per line,
    the codon_generator class can be initialized from a codon table with any number
    of entries per line. Additionally, this format allows for comment lines beginning with
    the '#' character, which will be ignored by the parser.

    NOTE: When first selecting a codon usage table from the Codon Usage Database, 
    one must select a genetic code to use in order to get values for 'amino_acid' 
    and 'fraction' to appear in the table.

    As this class was originally designed for generating random codons according to the
    biases of Sulfolobus solfataricus P2, if no codon usage table is provided
    when constructing a new codon_generator object, the codon usage table for 
    S. solfataricus P2 will be read from the Codon Usage Database, if an internet 
    connection is available.

    Each instance of a codon_generator object contains three attributes:
        self.codon_table    - A list of tuples (codon, amino_acid) where 'amino_acid'
                              is the amino acid coded for by 'codon'. This list contains
                              one entry for each possible codon.

        self.codon_cdf      - A list of floats, one for each codon, where 
                              codon_cdf[index] is the sum of the probabilities of each
                              codon in codon_table[0:index+1]. For example, codon_cdf[0]
                              is simply the probability of the codon in codon_table[0]
                              whereas codon_cdf[1] is the sum of the probabilities of the
                              codons in codon_table[0] and codon_table[1].

        self.aa_table       - A dict with single letter amino acid codes as keys, and
                              a pair of lists as values, one pair of lists per amino acid. 
                              That is, aa_table[aa] is the tuple (codons,aa_cdf) where
                              'codons' is a list of strings where each string is a codon
                              for the amino acid 'aa', and 'aa_cdf' is similar to 
                              'codon_cdf' defined above, but using codon probabilities
                              normalized by the total number of codons encoding 'aa'
                              that were observed, instead of normalized by the total
                              number of codons observed.

    The codon_generator class defines a constructor, two public methods, and a private
    method as follows:
        codon_generator     - An initializer method that optionally takes in an input 
                              stream containing a codon usage table in the format 
                              described above. This data is used to determine 
                              probabilities. If no codon usage table is provided, one 
                              for S. solfataricus P2 will be read from the web.

        get_random_codon    - Takes in a set of disallowed protein sequence characters,
                              and returns a random codon for an allowed residue according
                              to the provided codon usage table.

        get_random_aa_codon - Takes in a protein sequence character, and returns a 
                              random codon for that character, weighted by the codon
                              usage table.
       
        __check_cdfs        - A private method used internally to make sure that the 
                              calculated probabilities sum to 1.0, within a reasonable
                              rounding error.

    For more information on the methods provided, please see the docstring for the 
    individual function.
    """
    def __init__(self,input_stream=None):
        """Construct a new codon_generator object from a codon usage table.

        Inputs:
            input_stream - An input stream containing codon usage table data in the
                           format described in the codon_generator class docstring.
                           If not provided, a codon usage table for S. solfataricus P2
                           will be read from the web if possible.
        Outputs:
            codon_gen    - A codon_generator object initialized to codon probabilities
                           determined from the given codon usage table.

        This method is used to create a new codon_generator object using probabilities
        determined from a codon usage table either provided by the user or read from 
        the web.
        """
        #First we initialize the lists and dict that will become the fields for our new
        #codon_generator object, as well as 'codon_data' to temporarily store codon entries
        codon_table = []
        aa_table = dict()
        codon_data=[]
        codon_cdf = []
        total_codons=0.0
        if input_stream is None:
            #If no input stream provided, read from web and read lines until we get to 
            #the beginning of the codon table, marked by the HTML <PRE> tag on the line
            #right before the table starts.
            input_stream = urllib2.urlopen("http://www.kazusa.or.jp/codon/cgi-bin/"
                                   "showcodon.cgi?species=273057&aa=1&style=N")
            for line in input_stream:
                if line.startswith("<PRE>"): break
        for line in input_stream:
            #Want to skip over comment lines, or break if done reading the table from web.
            #As the 'number' field is often in parentheses, use the following re to split
            #the line into strings of non-whitespace, non-parentheses characters to get the
            #individual fields to iterate over. 'codon' will be used to keep track of which
            #codon we're currently reading data for.
            if line.startswith("</PRE>"): break
            if line.startswith('#'): continue
            elements = re.findall("[^\s()]+",line)
            codon=""
            for num,element in enumerate(elements):
                if num % 5 == 0:
                    if codon != "":
                        #Each entry has 5 fields, hence store a new entry every five 
                        #fields, except the first time through, don't want to store
                        #an empty codon. We pull out the individual fields from the
                        #temporary data list 'codon_data', and store the approriate
                        #values. For now, cdf lists contain cumulative sums of 
                        #raw counts of occurences and will be normalized after all 
                        #the data is read.
                        codon,aa,fraction,per_thousand,number = codon_data
                        number = int(number)
                        total_codons += number
                        codon_table.append((codon,aa))
                        codon_cdf.append(total_codons)
                        if aa not in aa_table:
                            #If no entry for amino acid 'aa' exists yet, add a tuple of 
                            #lists with the right values for the current amino acid.
                            aa_table[aa] = ([codon],[number])
                        else:
                            #If an entry already present, we want to update the total codon
                            #count for this amino acid, and store the codon and new total.
                            last_aa_total = aa_table[aa][1][-1]
                            aa_total = last_aa_total + number
                            aa_table[aa][0].append(codon)
                            aa_table[aa][1].append(aa_total)
                    #Update 'codon' to the new codon entry we're about to begin, and reset 
                    #the temporary data container 'codon_data'
                    codon=element
                    codon_data=[]
                #At each step, we just add the field value to our temporary data container
                #storing the data until we reach the beginning of a new entry.
                codon_data.append(element)
            #Need to add entries at the end of a line as well, hence the duplication of code below
            if codon != "":
                codon,aa,fraction,per_thousand,number = codon_data
                number = int(number)
                total_codons += number
                codon_table.append((codon,aa))
                codon_cdf.append(total_codons)
                if aa not in aa_table:
                    aa_table[aa] = ([codon],[number])
                else:
                    last_aa_total = aa_table[aa][1][-1]
                    aa_total = last_aa_total + number
                    aa_table[aa][0].append(codon)
                    aa_table[aa][1].append(aa_total)
        #After we've finished reading in all the data, just need to normalize raw counts
        #to get probabilities, either normalizing by the total number of codons observed
        #in the case of 'codon_cdf' or by the total number of codons for a specific amino
        #acid observed in the case of each 'aa_cdf'.
        for index in xrange(len(codon_cdf)):
            codon_cdf[index] /= total_codons
        for key in aa_table:
            aa_codons,aa_cdf = aa_table[key]
            aa_total = float(aa_cdf[-1])
            for index in xrange(len(aa_cdf)):
                aa_cdf[index] /= float(aa_total)
        #Make sure all probabilities sum to 1 within rounding error, and store data in
        #appropriate fields.
        self.codon_table = codon_table
        self.codon_cdf = codon_cdf
        self.aa_table = aa_table
        self.__check_cdfs()



    def __check_cdfs(self):
        """Checks to make sure probabilities sum to 1.0 within rounding error.

        __check_cdfs() takes no inputs and has no output, it simply makes sure that
        the probabilities described in 'self.codon_cdf', and in the 'aa_cdf's stored
        in 'self.aa_table', sum to 1.0 within a reasonable rounding error, defined as
        1E-10.
        
        If any probabilities sum to a value outside of the acceptable range, an error 
        will be printed and this function will call sys.exit(1). If the sum is within 
        the tolerance range, the final element of the cdf lists will be set to 1.0 to 
        prevent the possbibility of random.random() producing a number that causes 
        bisect.bisect(cdf_list, random_num) to return an index off the end of the cdf_list.
        """
        #For both codon_cdf and aa_cdf, the last element gives the sum of all 
        #probabilities, hence we just check self.codon_cdf[-1] and aa_cdf[-1]. If the
        #sum is within the tolerance range, we set the final element of the cdf list to
        #1.0.
        tolerance = 0.0000000001
        final_cdf = self.codon_cdf[-1]
        if abs(1.0 - final_cdf) > tolerance: 
            print("ERROR: Codon probabilities do not sum to 1.0", file = sys.stderr)
            sys.exit(1)
        else:
            self.codon_cdf[-1] = 1.0
        for aa_codons,aa_cdf in self.aa_table.itervalues():
            #For each amino acid, make sure codon probabilities are acceptable.
            final_cdf = aa_cdf[-1]
            if abs(1.0 - final_cdf) > tolerance:
                print("ERROR: Amino acid codon probabilities do not sum to 1.0 for "
                      "an amino acid." , file = sys.stderr)
                sys.exit(1)
            else:
                aa_cdf[-1] = 1.0


                     

    def get_random_codon(self,not_allowed=[]):
        """Return a random codon according to the biases in the codon usage table.

        Inputs:
            not_allowed - An optional iterable supporting the 'in' operator containing 
                          protein sequence characters for which codons are not desired.
                          If a generated codon encodes a character in 'not_allowed' it 
                          will be discarded and a new codon will be generated. A common 
                          use of this input option is to disallow stop codons by using:
                              codon_gen.get_random_codon(not_allowed="*")
                          Defaults to an empty set.
        Outputs:
            codon       - A random codon generated according to the codon usage table
                          used to initialize the codon_generator object calling this 
                          method. This codon will not code for a protein sequence 
                          character contained in 'not_allowed'.
  
        The get_random_codon function is used to generate a codon randomly according
        to the codon usage table used when the calling codon_generator object was 
        instantiated. If no iterable is provided for the not_allowed option, then
        any of the possible codons will be produced. 
        """
        #First we default to 'allowed_codon' being False, indicating we haven't generated
        #an allowed codon yet. Set to True when an allowed codon is generated, causing the
        #codon to be returned.
        allowed_codon = False
        while(not allowed_codon):
            #We use the library functions random.random() and bisect.bisect() to generate
            #a random number and translate it to an index in codon_table.
            rand = random()
            index = bisect(self.codon_cdf,rand)
            codon,aa = self.codon_table[index]
            if aa in not_allowed: continue
            else: allowed_codon = True
        return codon

    def get_random_aa_codon(self,aa):
        """Return a random codon, encoding the specified amino acid, according to the codon usage table.

        Input:
            aa    - Required argument giving a single letter amino acid code, or the 
                    stop symbol '*', specifying the amino acid to generate a random codon 
                    for.
        Output:
            codon - A codon encoding the amino acid 'aa' generated at random according to
                    the counts given in the codon usage table.

        The get_random_aa_codon function generates a random codon for the specified amino
        acid, 'aa'. The codon is randomly generated according to the counts given in the 
        codon uage table used when the calling codon_generator object was instantiated.
        """
        #Check to make sure the specified amino acid code was present in the codon table
        #and error out and exit if so.
        if aa not in self.aa_table:
            print("ERROR: Single-letter amino acid code provided not seen in codon "
                  "table: %s" % (aa), file=sys.stderr)
            sys.exit(1)
        #We use the library functions random.random() and bisect.bisect() to generate
        #a random number and translate it to an index in aa_codons, the list containing
        #the codons for the specified amino acids.
        rand = random()
        aa_codons,aa_cdf = self.aa_table[aa]
        index = bisect(aa_cdf,rand)
        codon = aa_codons[index]
        return codon
        

