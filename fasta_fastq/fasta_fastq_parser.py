#! /usr/bin/env python3.3
"""This module contains generators which yield FASTA or FASTQ sequences one at a time from a given input stream.

fasta_fastq_parser.py
BME 205 Fall 2014, Programming Assignment #2
October 17th, 2014
Robert Calef (rcalef@ucsc.edu and calef@soe.ucsc.edu)

This module defines two classes, 'fasta_seq' for representing FASTA sequences, and 'fastq_seq' for representing 
FASTQ sequences, and are described below:

fasta_seq: This class has no functions, and 3 fields:
   self.sequence - contains the FASTA sequence itself as a string
   self.identifier - the ID of the specific sequence as a string
   self.comment - contains any comments from the sequence ID line as a string

fastq_seq: This class also has no functions, and contains all 3 fields described above with the same
           naming convention, along with an additional 'scores' field to store quality scores:
   self.sequence - contains the FASTA sequence itself as a string
   self.identifier - the ID of the specific sequence as a string
   self.comment - contains any comments from the sequence ID line as a string
   self.scores - a list of integers, one quality score for each character in self.sequence

The above classes are used 
"""
from __future__ import print_function
import sys
import re
import operator
import argparse

#The below globally defined lists serve as default alphabets for the respective parsers
fastq_allowed_bases=['A','C','G','T','N','U','K','S','Y','M','W','R','B','D','H','V','-']
fasta_allowed_bases=['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','U','V','W','X','Y','Z','-','*']

class fastq_seq:
    """A convenience class for storing a FASTQ sequence and its associated data.

    fastq_seq is a simple convenience class containing no methods gathering 
    together the variables used to store a FASTQ sequence, quality scores, 
    sequence ID, and any other comments. The fields contained in each instance 
    of the class are described below:
   
        self.sequence   - A string containing the FASTQ sequence itself, with each 
                          symbol of the original sequence as a character.

        self.identifier - A string containing the sequence ID that was on the header
                         line of the original FASTQ entry. The sequence ID is 
                         defined as everything on the header line after the beginning 
                         '@' and up to the first whitespace character or comma. A
                          sequence with no ID will have an empty string in this field.

        self.comment    - A string containg the comment from the header line of the
                          original FASTA entry. The comment is defined as everything
                          on the header line after the first whitespace character or
                          comma. An entry with no comment will have an empty string 
                          as a comment.

        self.scores    -  A list of integers representing the Phred quality scores 
                          for each symbol in the sequence. Specifically, there is 
                          guaranteed to be one score for each character in 
                          self.sequence, where self.scores[i] is the quality score for 
                          self.sequence[i].
    """
    def __init__(self,comment,identifier,sequence,qual_scores):
        self.comment = comment
        self.identifier = identifier
        self.sequence = sequence
        self.scores=qual_scores

class fasta_seq:
    """A convenience class for storing a sequence in FASTA format.

    fasta_seq serves as a simple convenience class to hold a FASTA sequence and all
    of its associated data. This class contains no functions. The fields of each
    instance of the class are described below:

        self.sequence   - A string containing the FASTA sequence itself, each symbol
                          in the original sequence as a character.

        self.identifier - A string containing the sequence ID that was on the header
                          line of the original FASTA entry. The sequence ID is
                          defined as everything on the header line after the beginning
                          '>' and up to the first whitespace character or comma. A
                          sequence with no ID will have an empty string in this field.

        self.comment    - A string containg the comment from the header line of the
                          original FASTA entry. The comment is defined as everything
                          on the header line after the first whitespace character or
                          comma. An entry with no comment will have an empty string
                          as a comment.

    """
    def __init__(self,comment,identifier,sequence):
        self.comment = comment
        self.identifier = identifier
        self.sequence = sequence


def header_split(header):
    """A utility function taking in a FASTA/FASTQ header line and returning the ID and comment.

    Input:
        header  - a FASTA or FASTQ header line as a string

    Output:
        (seq_ID,comment) - A tuple containing the sequence ID as its first element
                           and any associated comments as its second element.

    header_split strips off any trailing newline, as well as the starting '>' or '@' 
    characters, and then finds the first occurence of a whitespace character or a comma 
    in the input text. The text is then split around this separator, discarding the 
    separator and storing everything before the separator as 'seq_ID' and everything
    after the separator as 'comment', with both of these variables returned as a tuple.
    """
#Remove any trailing white space or new lines
    header = header.strip()
#Use the regular expression "[\s,]" which matches any whitespace character or a comma 
#to get index of first whitespace or comma.
    match = re.search("[\s,]",header)
#If no match occurs, return the header without the initial character as the sequence ID
    if match is None: return(header[1:],"")
#Otherwise split the header about the beginning of the first match to get ID and comment
    seq_id=header[1:match.start(0)]
    comment=header[match.start(0)+1:]
    return (seq_id,comment)

def read_fastq(input, phred, alphabet=None,filter_seqs=False,ignore_case=True):
    """A generator taking in an input stream and Phred score offset, and yielding FASTQ sequences.

    Inputs:
      Required:
        input - An input stream, such as an opened file or sys.stdin, containing data in FASTQ
                format. Lines will be read from this input to generate fastq_seq objects.

        phred - An integer specifying the Phred score offset for the quality scores in the 
                FASTQ data, typically 33 or 64 for Phred33 and Phred64 formats respectively.
      Optional:
        alphabet    - A set of characters in any iterable supporting the "in" operator defining
                      the alphabet of allowable characters in the FASTQ sequence itself. If no
                      alphabet is specified, the default FASTQ alphabet defined at the beginning
                      of the file will be used, containing all standard nucleotide and degenerate
                      nucleotide symbols, as well as '-' as the gap character. Any  characters not 
                      in the specified alphabet will be discared from the sequence. 
                      WARNING: this will not remove any quality information for discarced
                      characters, which could lead to different lengths of the sequence and 
                      the quality score string, causing the sequence to not be yielded.

        filter_seqs - A boolean specifying whether or not the parser should act as a filter as
                      well. Functionally, all this means is that if filter_seqs is False,  
                      warnings will be printed to stderr for each non-whitespace character
                      encountered that isn't in the defined alphabet. If filter_seqs is true
                      then these errors will be suppressed. Be very careful using filter_seqs=True 
                      with a user-defined alphabet.

        ignore_case - A boolean specifying whether or not to ignore the case of letters in the
                      FASTQ sequence when checking if in the defined alphabet, defaults to True.
                      Note that this only affects membership in the alphabet, the sequence itself
                      will be yielded with all cases preserved as in the original sequence.
    
    Output:
        fastq_seq   - Each FASTQ entry in the input stream will be yielded one at a time as a 
                      fastq_seq object as defined above. Any empty IDs, comments, sequences or
                      quality score strings will be stored as empty strings or an empty list
                      in their respective fields.

    read_fastq will parse an input stream containing data in FASTQ format, and yield each FASTQ entry
    as a fastq_seq object, one entry at a time, and terminating upon reaching the end of file. Any
    FASTQ entries with fewer quality symbols than bases will produce buggy output, but FASTQ entries
    with more quality symbols than bases will print an error, and truncate the quality score data to
    match the length of the sequence data.
    """
#bases_read and scores_read used in the parsing loop to keep track of the number of bases and quality
#scores read for each sequence respectively, and report an error if not equal.
    bases_read=scores_read=0
#in_scores and in_sequence are two booleans used in the parsing loop to define where in the file
#we currently are, whether we're reading in sequences, scores, or are between FASTQ entries.
    in_scores=in_sequence=False
#Initialize variables to hold sequence ID, comment, scores, and actual sequence as we parse each FASTQ
#entry to build a fastq_seq object.
    seq_id=comment=sequence=""
    scores=[]
#Check for user-defined alphabet, if none, then use default defined at the top of the file.
    if alphabet is None:
        alphabet = fastq_allowed_bases
    for line in input:
#If we've reached a new line while reading scores, and have read a score for each base we saw in
#the actual sequence, then we're done with this FASTQ entry (assuming the quality score string
#is actually the right length, this will be buggy if too few quality scores), so yield the 
#fastq_seq object and set in_sequence=in_scores=False to indicate we're between entries.
        if in_scores and bases_read==scores_read:
            yield fastq_seq(comment,seq_id,sequence,scores)
            in_sequence=False
            in_scores=False
#If not reading sequence, then we're not in an entry, so check if this line is a header of a new
#entry, and if not, keep going until we find the beginning of the next entry.
        if not in_sequence:
            if not line.startswith('@'):
                continue
#If actually a header line, then get sequence ID and comment, reset temporary sequence and scores
#containers, as well as bases and scores counters, and indicate that we're now in a sequence.
            id_comment = header_split(line)
            seq_id = id_comment[0]
            comment = id_comment[1]
            sequence = ""
            scores = []
            in_sequence = True
            bases_read = scores_read = 0
            continue
        if line.isspace(): continue
#If we run in to a line starting witha '+' while not reading scores (also while reading sequence as
#if we're not in_sequence, we'll continue in the above conditional block until we are) then this
#signals the start of score information on the next line.
        if not in_scores and line.startswith('+'):
            in_scores = True
            continue
#Above code largely for logic in terms of where in a FASTQ entry we are, now we read each character
        for char in line:
#If white space, discard
            if char.isspace(): continue
#Like above, at this point we're guaranteed to be reading sequences or scores, so if not scores,
#then check if sequence symbol is in alphabet, checking ignore_case for appropriate action,
#and print warning message if not found, only discard non-whitespace symbols not in alphabet
#if filter_seqs is set to True by user.
            if not in_scores:
                if ignore_case: check = char.upper()
                else: check = char
                if check not in alphabet:
                        if not filter_seqs:
                            print("WARNING Encountered sequence symbol not"
                            " in alphabet: %s\nSequence & position: %s | %d"
                            % (char,seq_id,pos),file=sys.stderr)
                        else:
                            continue
                sequence += char
                bases_read += 1
#If reading scores, check if score is non-negative, if so, add to scores list.
            else:
                score=ord(char)-phred
                if score < 0: continue
#If scores_read ==bases_read at this point, then we've read more scores than bases and discard
#any excess quality scores.
                if scores_read == bases_read:
                    print("ERROR Number of bases and number of quality "
                          "characters not the same: %s\n Truncating "
                          "quality scores to match number of bases."
                          % (seq_id),file=sys.stderr)
                    break
                scores.append(score)
                scores_read +=1
#If we reach the end of file while reading an entry, then return last entry
#(in_scores implies in_sequence by design).
    if in_sequence:
        yield fastq_seq(comment,seq_id,sequence,scores)

def read_fasta(input,alphabet=None,filter_seqs=False,case_sensitive=False):
    seq_id = ""
    comment = ""
    sequence=""
    in_sequence = False
    if alphabet is None: 
        alphabet = fastq_allowed_bases
#Now iterate over input until a new sequence is found
    for line in input:
#If first character of line is '>', then we've encountered the
# start of a new sequence entry, need to set this line to be
#curr_line and return the current sequence.
        if in_sequence and line.startswith('>'):
            yield fasta_seq(comment,seq_id,sequence)
            in_sequence = False
        if not in_sequence:
            if not line.startswith('>'): continue
            id_comment = header_split(line)
            seq_id = id_comment[0]
            comment = id_comment[1]
            sequence=""
            in_sequence = True
            continue
        if line.isspace(): continue
#Check if valid character
        for char in line:
#If white space, continue
            if char.isspace(): continue
            if case_sensitive: check = char
            else: check = char.upper()
            if check not in alphabet:
                if filter_seqs:
                    continue
                else:
                    print("WARNING Encountered sequence symbol not"
                    " in alphabet: %s\nSequence & position: %s | %d"
                    % (char,seq_id,pos),file=sys.stderr)
            sequence += char
#If we reach the end of file, then return last entry
    if in_sequence:
        yield fasta_seq(comment,seq_id,sequence)

def read_fasta_with_quality(fasta_input, qual_input):
#First get FASTA reader and initialize variables used throughout loop
    fasta_seqs=read_fasta(fasta_input)
    in_score = False
    scores=[]
    for line in qual_input:
#If currently reading a score, and find a line beginning with '>'
#then get the next FASTA sequence, and yield along with scores
#as a new FASTQ sequence
        if in_score and line.startswith('>'):
#Need to make sure FASTA file contains same number of sequences as qual file
            try:
                fasta_seq = next(fasta_seqs)
            except StopIteration:
                print("Number of FASTA sequences not "
                    "equal to number of quality "
                    "sequences", file=sys.stderr)
                return
#Also need to make sure FASTA sequence and qual sequence are same length
#if not, then yield None to signify a bad sequence
            if len(scores) != len(fasta_seq.sequence):
                print("FASTA sequence: ", fasta_seq.identifier,"\n",
                      "is not same length as quality score sequence"
                      " with same identifier. Omitting sequence."
                      ,file=sys.stderr)
#Else yield the new sequence and reset in_score to indicate we need
#to start a new sequence.
            else: yield fastq_seq(fasta_seq.comment,fasta_seq.identifier,fasta_seq.sequence,scores)
            in_score=False
        if not in_score:
            if not line.startswith('>'): continue
            id_comment = header_split(line)
            qual_id = id_comment[0]
            scores=[]
            in_score=True
            continue
#If we're reading a block (in_score==True) and the line does not begin
# with '>', then split the line into quality scores, check each to
#make sure they're a valid digit, and append to current score sequence
        for num in line.strip().split():
            if(num.isdigit()):
                score=int(num)
                if score >= 0: scores.append(score)
#When we reach end of file, we still have to yield the final sequence,
#hence the duplication of code below
    if in_score:
        try:
            fasta_seq = next(fasta_seqs)
        except StopIteration:
            print("Number of FASTA sequences not "
                  "equal to number of quality "
                  "sequences", file=sys.stderr)
            return
        if len(scores) != len(fasta_seq.sequence):
            print("FASTA sequence: ", fasta_seq.identifier,"\n",
                  "is not same length as quality score sequence"
                  " with same identifier. Omitting sequence."
                  ,file=sys.stderr)
        else: yield fastq_seq(fasta_seq.comment,fasta_seq.identifier,fasta_seq.sequence,scores)
    else:
        try:
            fasta_seq = next(fasta_seqs)
        except StopIteration:
            pass
        else:
            print("Number of FASTA sequences not "
                  "equal to number of quality "
                  "sequences", file=sys.stderr)
            return

 

def main(args):
    for arg in args:
        print(arg)
    fasta_file = open(args[1],'r')
    qual_file = open(args[2],'r')
    for sequence in read_fasta_with_quality(fasta_file,qual_file):
        if sequence is None: continue
        sys.stdout.write("\n@%s   %s" % (sequence.identifier,sequence.comment))
        for i in range(len(sequence.sequence)):
            if i % 80 == 0: print
            sys.stdout.write(sequence.sequence[i])
        print("\n+")
        i=0
        for score in sequence.scores:
            if i % 80 == 0: print
            sys.stdout.write(chr(score+33))
            i+=1
        print
    print("\n")
#    print_dict(words, order)



if __name__ == "__main__":
    sys.exit(main(sys.argv))
