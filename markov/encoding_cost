#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse
import string
from math import log
from collections import Counter
from operator import itemgetter

from fasta_fastq_parser import read_fasta
from markov import generate_all_kmers, kmers_from_sequence, check_alphabet


def parse_arg():
    """Parses command line arguments, returning appropriate values in a Namespace object.

    parse_arg() takes no inputs, as the parse_args() method accesses the command line 
    arguments directly. This function constructs an ArgumentParser object with various
    options specified in the --help message of this program.

    The only required argument is the positional 'fasta' argument specifying the file
    from which to read sequences and calculate encoding costs. For a detailed description
    of each option, see the --help message.
    """
    #Construct initial ArgumentParser object
    argparser = argparse.ArgumentParser(description = __doc__)
    #Default alphabet is an empty string as we will add characters as we read counts
    argparser.add_argument('-a','--alphabet',action='store', 
        default="",
        help=("Specifies the alphabet of valid characters in the"
        " FASTA sequences, any other characters will be discarded. Be careful "
        "when specifying an alphabet, as dropping characters from a sequence can "
        "artificially introduce new kmers to the sequence. By default, "
        "the alphabet will be detected from the list of counts read from stdin, "
        "any character occuring in the list of counts will be added to the "
        "alphabet. WARNING: this program does not have the functionality to "
        "detect start and stop characters, if the counts contain a start character "
        "other than '^' or a stop character other than '$', they must be specified "
        "using the --start and --stop options. If an alphabet is specified, "
        "and non-alphabet characters are detected in the input counts, then an "
        "error will be printed and the program will exit."))
    argparser.add_argument('--start', action = 'store',default='^',
        help =("Specifies the character to be used to represent"
        " the beginning of a sequence. If not specified, '^' will be used."))
    argparser.add_argument('--stop', action = 'store',default='$',
        help =("Specifies the character to be used to represent"
        " the end of a sequence. If not specified, '$' will be used."))
    argparser.add_argument('-t','--train',action='store',help=("Optionally specify"
        " the file used to generate the counts to train the model. The program will "
        "then calculate and output encoding data for the training data as well."))
    argparser.add_argument('fasta', action='store',help=("A mandatory argument "
        "specifying the FASTA file containing data to calculate encoding costs for."))
    opts = argparser.parse_args()
    #Next we want to make sure the specified alphabet does not contain the start or stop
    #characters, and error if so. Also print warnings if non-printable or whitespace
    #characters are specified in the alphabet.
    check_alphabet(opts.alphabet,opts.start,opts.stop)
    return opts
    

def read_counts(start,stop,input_file,alphabet=""):
    """Reads tab-delimited kmer,count pairs, one per line, from an input stream.

    Inputs:
        start      - The character representing the start of a sequence.
        stop       - The character represeting the end of a sequence.
        input_file - The input stream from which to read count data.
        alphabet   - An optional set object containing the alphabet of valid characters
                     for count data. If any non-alphabet characters are found in the 
                     kmers, an error will be printed and the program will exit.
    Outputs:
        Returns a tuple (counts,k):
          counts   - A Counter object containing the parsed count data, given some kmer
                     "AAA" in the input, counts["AAA"] will be the count associated with
                     that kmer.
          k        - An integer giving the length of the kmers.

    read_counts is a short utility function used to read in the kmer-count data used to
    build a Markov model. Count data must be formatted as tab-delimited kmer,count pairs,
    one pair per line, for example:

        AAA	2
        AAB	4
        ...

    although the data need not be sorted. Count data will be returned in a Counter object
    in a tuple along with k, the length of the kmers in the parsed data.
    """
    #First initialize the Counter object used to store the parsed counts
    counts = Counter()
    #specified_alphabet is used to alter non-alphabet character handling depending on
    #whether or not the user specified an alphabet,default alphabet is the empty string
    specified_alphabet=False
    if alphabet not == "": specified_alphabet=True
    for line in input_file:
        #First we split the line into the kmer and the count, and then check each
        #kmer's letters for membership in the alphabet
        kmer_count=line.split()
        kmer=kmer_count[0]
        count=kmer_count[1]
        for letter in kmer:
            if letter not in alphabet:
                #Error out if non-alphabet character with user-specified alphabet
                if specified_alphabet:
                    print("Counts data contains a character not in the specified "
                          "alphabet: %s" % (letter),
                           file=sys.stderr)
                    sys.exit(1)
                else: alphabet.add(letter)
        counts[kmer] = int(count)
    #Store the length of the kmers, and get rid of the start and stop characters in
    #the alphabet, as the get put in the alphabet when reading counts that should
    #contain kmers with start and stop characters. We don't want start and stop
    #characters in the alphabet, as these are not valid characters in raw sequences.
    k = len(kmer)
    alphabet.discard(start)
    alphabet.discard(stop)
    return (counts,k)


def make_conditional_log_probs(counts):
    """Constructs the conditional log probability table needed for an order N > 0 Markov model.

    Inputs:
        counts    - A dict-like object containing (kmers,count of occureces) as (key,value) 
                    pairs. These counts are assumed to already contain pseudocounts, and 
                    will be used directly to construct the probability table.
    Outputs:
        log_probs - A dict object containing conditional log probabilities for an order
                    N Markov model. The order of the model is k-1, where k is the kmer 
                    length. Given some kmer, the entry in this table will contain the
                    probability of seeing the k-th character given the k-1 characters
                    preceding it. E.g. log_probs["ABC"] is the base 2 logarithm of the
                    conditional probability of seeing 'C' given that the preceding two
                    letters were 'AB'.

    make_conditional_log_probs is a simple transformation function that takes in a dict
    of kmer-count pairs, and outputs a dict containing kmers as keys and the base 2 
    logarithm of the conditional probability of seeig the k-th character of the kmer
    given the k-1 characters preceding it. 

    Probabilities are constructed by simply dividng each kmer's count of occurences by 
    the sum of all kmer occurences with the same k-1 beginning characters. For example,
     if constructing a log-prob table for an order 2 Markov model, the conditional 
    probability of seeing an "A" given seeig "AA" right before is represented by the 
    3-mer "AAA". This conditional probability would be calculating by dividing the 
    number of occurrences of "AAA" by the sum of all 3-mer occurences beginning with the
    2-mer "AA", that is the sum of all conditional probabilities of seeing any character
    following "AA", including the stop character, which must sum to 1.
    """
    #We sort our list of counts (assumed to already contain pseudocounts) so as to group
    #kmers beginning with the same (k-1)-mer together. 
    kmers = sorted(counts)
    #context will be used to store the current (k-1)-mer, or context, of the kmers being
    #processed. This is used to detect when we've reached the end of a context and have 
    #to go back and normalize the appropriate kmer counts.
    context=None
    #Initialize a empty dict object to store the log probabilities to be calculated
    log_probs = dict()
    #num_kmers to keep track of the number of kmers in a given context, used to 
    #backtrack over kmers in a give context.
    num_kmers=0
    for itor,kmer in enumerate(kmers):
        #Get the context of the current kmer, if it is not the same as the context
        #seen in previous kmer, then we need to go back, normalize, and repeat.
        curr_context=kmer[0:-1]
        #Note that context is initialized to None, so the following conditional is entered
        #on the first iteration of this for loop
        if context != curr_context:
            #Starting a new context, store the new context
            context=curr_context
            #For each kmer in the previous context, normalize the count
            for context_kmer in kmers[itor-num_kmers:itor]:
                #conditional probability is the count of the kmer over counts of all 
                #kmers in that context
                prob = float(counts[context_kmer])/float(context_sum)
                #After we get the probability, store the base 2 log in the log prob table
                log_probs[context_kmer]= -log(prob,2)
            #Finally, reset the context sum and number of kmers, and begin the 
            #next context
            context_sum=0
            num_kmers=0
        num_kmers += 1
        context_sum += counts[kmer]
    #After exiting the preceding for loop, we still have to normalize the final context
    #hence the duplication of code below.
    itor += 1
    for context_kmer in kmers[itor-num_kmers:itor]:
        prob = float(counts[context_kmer])/float(context_sum)
        log_probs[context_kmer]= -log(prob,2)
    return log_probs



def make_log_frequency_table(counts):
    """Constructs the character probability table used for an order 0 Markov model.

    Inputs:
        counts    - A dict-like object containing (kmers,count of occureces) as (key,value)
                    pairs. These counts are assumed to already contain pseudocounts, and
                    will be used directly to construct the probability table.
    Output:
        log_probs - A dict object containing the base 2 logarithm of individual 
                    character probabilities for an 
                    order 0 Markov model, or a model looking at character probability
                    independent of any preceding characters. Each (key,value) pair in
                    the dict is a (character, log_2(probability)) pair where the probability
                    for a character is the count of that character's occurences divided
                    by the total number of characters.

    make_log_frequency_table transforms a dict of (character, count of occurences) 
    key-value pair to a dict of (character, log_2(probability of occurence)) key-value 
    pairs. The probability for each character is simply it's number of occurences over 
    the total number of characters.
    """
    #First we get the total number of characters seen.
    total_counts=0
    #counts.items() returns an iterator over the (key,value) tuples in the table.
    for kmer_counts in counts.items():
        total_counts += kmer_counts[1]
    #Then we initialize the emoty log-prob table, and populate it.
    log_probs = dict()
    for kmer_counts in counts.items():
        prob = float(kmer_counts[1])/float(total_counts)
        log_probs[kmer_counts[0]] = -log(prob,2)
    return log_probs


def add_pseudocounts(pseudo,counts,alphabet,start,stop,k):
    """A simple utility function to add uniform pseudocounts for all possible kmers to a Counter object.

    Inputs:
        pseudo   - The pseudocount value to be added to all possible kmers.
        counts   - A Counter object containing (kmer,count of occurences) as (key,value) 
                 pairs, counts in this object will be incremented directly.
        alphabet - The alphabet of valid sequence characters, not including the start
                   or stop characters.
        start    - The character used to represent the beginning of a sequence.
        stop     - The character used to represent the end of a sequence.
        k        - The length of kmers for which pseudocounts are being generated.
    Outputs:
        none, simply increments the values in 'counts', and adds new entries for 
        unobserved kmers
    
    add_pseudocounts is a simple function that leverages generate_all_kmers to increment 
    counts for all possible kmers in a Counter object 'counts' by some set value 'pseudo'.
    Entries will be added for unobserved kmers that are valid, as defined in 
    generate_all_kmers.
    """
    for kmer in generate_all_kmers(k,alphabet,start,stop):
        counts[kmer] += 1

def get_sequence_coding_cost(sequence,alphabet,log_probs,start,stop,k):
    """

    Inputs:
        pseudo   - The pseudocount value to be added to all possible kmers.
        counts   - A Counter object containing (kmer,count of occurences) as (key,value)
                 pairs, counts in this object will be incremented directly.
        alphabet - The alphabet of valid sequence characters, not including the start
                   or stop characters.
        start    - The character used to represent the beginning of a sequence.
        stop     - The character used to represent the end of a sequence.
        k        - The length of kmers for which pseudocounts are being generated.


    """
    cost=0.0
    if k==1:
        sequence = string.upper(sequence) + stop
    else:
        sequence= (start * (k-1)) + string.upper(sequence) + stop
    alphabet.add(start)
    alphabet.add(stop)
    for kmer in kmers_from_sequence(sequence,k):
        cost += log_probs[kmer]
    return cost

def get_file_coding_cost(fasta_file,log_probs,alphabet,start,stop,k):
    cost = 0.0
    tot_chars = 0
    num_sequences = 0
    for fasta_seq in read_fasta(fasta_file,alphabet,filter_seqs=True,ignore_case=True):
        seq_data = get_sequence_coding_cost(fasta_seq.sequence,alphabet,log_probs,
                                             start,stop,k)
        tot_chars += len(fasta_seq.sequence)
        cost += seq_data
        num_sequences += 1
    return(cost,num_sequences,tot_chars)

def print_coding_cost_data(file_name,total_cost,num_sequences,num_chars,k):
    print("Order %d Markov model\nTest file: %s\n"
          "Total encoding cost of test file: %g\nAverage bits/sequence: %g\n"
          "Avg bits/char: %g\n"
          % (k-1,file_name,total_cost,(total_cost/num_sequences),
          (total_cost/num_chars)))

def main(args):
    opts=parse_arg()
    alphabet=set(opts.alphabet)
    fasta_file = open(opts.fasta,'r')
    if opts.train is not None:
        train_file = open(opts.train,'r')
    (counts,k)=read_counts(alphabet,sys.stdin)
    alphabet.discard(opts.start)
    alphabet.discard(opts.stop)
    add_pseudocounts(k,counts,alphabet,opts.start,opts.stop)
    if k == 1:
        counts.pop(opts.start, None)
        log_probs = make_log_frequency_table(counts)
    else:
        log_probs = make_conditional_log_probs(counts,len(alphabet),k)
    (total_cost,num_sequences,num_chars)=get_file_coding_cost(fasta_file,log_probs,
                                              alphabet,opts.start,opts.stop,k)
    print_coding_cost_data(fasta_file.name,total_cost,num_sequences,num_chars,k)
    fasta_file.close()
    if opts.train is not None:
        (total_cost,num_sequences,num_chars)=get_file_coding_cost(train_file,
                                    log_probs,alphabet,opts.start,opts.stop,k)
        print_coding_cost_data(train_file.name,total_cost,num_sequences,num_chars,k)
        train_file.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
