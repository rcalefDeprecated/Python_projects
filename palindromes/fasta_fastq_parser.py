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
import string

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

def read_fastq(input_file, phred, alphabet=None,filter_seqs=False,ignore_case=True):
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
                      in the specified alphabet will be discarded from the sequence. 
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
                      in their respective fields, which are briefly described below for completness
                      of this function description:
                         fastq_seq.identifer - String containing sequence ID
                         fastq_seq.comment   - String containing any comments from the header line
                         fastq_seq.sequence  - String containing the sequence itself
                         fastq_seq.scores    - List of integers specifying quality scores

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
    for line in input_file:
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
                            " in alphabet: %s\nSequence: %s "
                            % (char,seq_id),file=sys.stderr)
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

def read_fasta(input_file,alphabet=string.ascii_uppercase,filter_seqs=False,ignore_case=True):
    """A generator taking in an input stream of FASTA formatted data, and yielding FASTA sequences.

    Inputs:
      Required:
        input - An input stream, such as an opened file or sys.stdin, containing data in FASTA
                format. Lines will be read from this input to generate fasta_seq objects.
      Optional:
        alphabet    - A set of characters in any iterable supporting the "in" operator defining
                      the alphabet of allowable characters in the FASTA sequence itself. If no
                      alphabet is specified, the default FASTA alphabet defined at the beginning
                      of the file will be used, containing all standard nucleotide and degenerate
                      nucleotide symbols, one letter amino acid codes, as well as '-' and '*' as 
                      the gap character and translation stop symbol respectively. Any  characters 
                      not in the specified alphabet will be discarded from the sequence. 

        filter_seqs - A boolean specifying whether or not the parser should act as a filter as
                      well. Functionally, all this means is that if filter_seqs is False,  
                      warnings will be printed to stderr for each non-whitespace character
                      encountered that isn't in the defined alphabet. If filter_seqs is True
                      then these errors will be suppressed. Be very careful using filter_seqs=True 
                      with a user-defined alphabet.

        ignore_case - A boolean specifying whether or not to ignore the case of letters in the
                      FASTA sequence when checking if in the defined alphabet, defaults to True.
                      Note that this only affects membership in the alphabet, the sequence itself
                      will be yielded with all cases preserved as in the original sequence.
    
    Output:
        fasta_seq   - Each FASTA entry in the input stream will be yielded one at a time as a 
                      fasta_seq object as defined at the beginning of the file. Any empty IDs, 
                      comments, or sequence will be stored as empty strings in their respective fields
                      which are briefly described below for completeness of this function description:
                         fasta_seq.identifer - String containing sequence ID
                         fasta_seq.comment   - String containing any comments from the header line
                         fasta_seq.sequence  - String containing the FASTA sequence itself

    read_fasta will parse an input stream containing data in FASTA format, and yield each FASTA entry
    as a fasta_seq object, one entry at a time, and terminating upon reaching the end of file. For
    the purposes of this parser, a FASTA entry is defined as:

        #header line
        >sequence_ID comments come after sequence ID
        #sequence data, possibly over multiple lines
        ACGTTACG...
        .
        .
        .

    where the sequence ID is defined as everything on the header line up to the first whitespace
    character or comma.
    """
#Initialize temporary container variables seq_id, comment, and sequence
#to store sequence ID, comments, and the sequence for each FASTA entry 
#as its built up line by line.
    seq_id = ""
    comment = ""
    sequence=""
#in_sequence is a boolean used to specify when the parser is currently
#reading lines thought to be sequence data, helps control the logic
#of the parsing loop.
    in_sequence = False
#If alphabet is None, then no user-defined alphabet -> use default.
    if alphabet is None: 
        alphabet = fastq_allowed_bases
    for line in input_file:
#If we reach a line beginning with '>' while reading sequences, then we've
#just hit the header for the next FASTA entry, so yield the current FASTA
#entry, and set in_sequence to False to indicate a new FASTA entry beginning.        
        if in_sequence and line.startswith('>'):
            yield fasta_seq(comment,seq_id,sequence)
            in_sequence = False
#If we're not parsing sequences, we want to continue until we reach the next
#header line, and once we reach a header line (starting with '>'), parse it
#to get the sequence ID and commment of the new FASTA entry, then finally
#indicate that we should now be reading in sequence data until another header
#line is encountered.
        if not in_sequence:
            if not line.startswith('>'): continue
            id_comment = header_split(line)
            seq_id = id_comment[0]
            comment = id_comment[1]
            sequence=""
            in_sequence = True
            continue
#Discard whitespace lines
        if line.isspace(): continue
#At this point in the loop, we are guaranteed to be reading sequence data, as
#if in_sequence is False, then we'll have continued in the previous conditional
#until reaching a new FASTA entry. So, start parsing characters.
        for char in line:
#If whitespace character, continue
            if char.isspace(): continue
#Check whether we're ignoring case, and then check if the appropriate case
#of the current letter is in the defined alphabet.
            if ignore_case: check = char.upper()
            else: check = char
            if check not in alphabet:
#If acting as a filter, do not print warning for characters not in defined alphabet
                if filter_seqs:
                    continue
                else:
                    print("WARNING Encountered sequence symbol not"
                    " in alphabet: %s\nSequence: %s"
                    % (char,seq_id),file=sys.stderr)
#If character passed all of the above checks, add to the sequence
            sequence += char
#If we reach the end of file while still building a new FASTA entry, then return it
    if in_sequence:
        yield fasta_seq(comment,seq_id,sequence)

def read_fasta_with_quality(fasta_input, qual_input,alphabet=None,filter_seqs=False,ignore_case=True):
    """A generator taking in input streams for FASTA and quality score data, and yielding FASTQ sequences.

    Inputs:
      Required:
        fasta_input - An input stream, such as an opened file or sys.stdin, containing data in FASTA
                      format. Lines will be read from this input to generate fasta_seq objects using
                      the read_fasta parser defined above.

        qual_input  - An input stream, such as an opened file or sys.stdin, containing quality score
                      data in a FASTA-like format described in detail below. Briefly, the quality
                      score file should have one entry per FASTA entry, in the same order as the
                      entries in the FASTA file, with whitespace separated integers specifying
                      quality scores on the line after the header line.
    
    Output:
        fastq_seq   - For each entry in the quality score input stream, a FASTA entry will be read from 
                      the FASTA input stream, and these will then be merged into a single fastq_seq  
                      object as defined at the beginning of the file. Any empty IDs, comments, sequence, 
                      or quality score string will be stored as empty strings or an empty list in their 
                      respective fields which are briefly described below for completeness of this 
                      function description:
                         fastq_seq.identifer - String containing sequence ID
                         fastq_seq.comment   - String containing any comments from the header line
                         fastq_seq.sequence  - String containing the sequence itself
                         fastq_seq.scores    - List of integers specifying quality scores

    read_fasta_with_quality takes in two input streams as its inputs, the first containing data in FASTA 
    format, and the second containing quality score data with the following specifications:

      -Each entry begins with a header line, exactly like a FASTA entry header, starting with '>'
       and containing the sequence ID and commment.
      -The following line (or lines) consist of whitespace separated integers, the i'th integer
       specifying the Phred quality score of the i'th symbol in the corresponding FASTA sequence
      -The quality score entries appear in the exact same order as the entries in the FASTA file
       i.e. Entry 1 in the quality file is the quality data for Entry 1 in the FASTA file, and so on

    An example entry:
  
       >sequence_ID comment, both same as corresponding entry in FASTA filei
       30 28 54 34 12 64 23 7 5 2 19 ...
       6 59 32 ...
       .
       .
       .
    Errors will be printed if there are either fewer or more FASTA entries than quality score entries,
    and pairs of FASTA and quality score entries that have mismatched sequence and quality score length
    will be discarded after an error is reported.  
    """
#First get FASTA parser and initialize scores, used to temporarily hold quality score data as
#its collected for a new fastq_seq object.
    fasta_seqs=read_fasta(fasta_input,alphabet,filter_seqs,ignore_case)
    scores=[]
#in_score indicates whether or not the parser is currently reading score data, or is between
#quality score entries
    in_score = False
    for line in qual_input:
#If currently reading a score, and find a line beginning with '>'then we've reached the end 
#of the current quality score entry, so we get ready to yield a fastq_seq object
        if in_score and line.startswith('>'):
#Need to make sure FASTA file contains same number of sequences as qual file, hence we try
#to get the next FASTA entry, if StopIteration is raised, then there were no more FASTA
#entries, and hence the number of FASTA entries was not the same as the number of quality
#score data entries.
            try:
                fasta_seq = next(fasta_seqs)
            except StopIteration:
                print("Number of FASTA sequences not "
                    "equal to number of quality "
                    "sequences", file=sys.stderr)
                return
            if not qual_id == fasta_seq.identifier:
                print("Mismatched FASTA entry and quality score entry:\n"
                      "   FASTA ID: %s\n"
                      "   Quality score ID: %s\n" % (fasta_seq.identifier,qual_id)
                      ,file=sys.stderr)
#We also need to make sure the FASTA sequence and the quality score sequence are the same length
#if not, then report an error and discard the FASTA/quality score pair.
            elif len(scores) != len(fasta_seq.sequence):
                print("FASTA sequence: ", fasta_seq.identifier,"\n",
                      "is not same length as quality score sequence"
                      " with same identifier. Omitting sequence."
                      ,file=sys.stderr)
#If all checks are passed, yield a new fastq_seq object constructed from our FASTA sequence data
#and the quality score data from the quality file. Set in_score to False to indicate that we're
#now looking for the next quality score entry.
            else: yield fastq_seq(fasta_seq.comment,fasta_seq.identifier,fasta_seq.sequence,scores)
            in_score=False
#If we're not currently reading score data, then we need to get the beginning of the next entry
#in the quality score file, hence continue until we reach a line starting with '>', and then
#parse sequence ID for error checking (make sure FASTA entry and quality entry have same ID).
        if not in_score:
            if not line.startswith('>'): continue
            id_comment = header_split(line)
            qual_id = id_comment[0]
            scores=[]
            in_score=True
            continue
#At this point in the loop, in_score is guaranteed to be True, as we continue in the above
#conditional until we reach the beginning of a new quality score entry, so start trying to
#parse quality scores. First strip trailing whitespace, and split the line about whitespace.
        for num in line.strip().split():
#Check if the non-whitespace string is a number, if so, convert to int and add to scores if
#the score is greater than 0 (Phred scores can't be negative).
            if(num.isdigit()):
                score=int(num)
                if score >= 0: scores.append(score)
#If we reach end of file while reading scores, then we still need to yield the final fastq_seq
#object, so try to get a final FASTA sequence and yield one more fastq_seq object.
    if in_score:
        try:
            fasta_seq = next(fasta_seqs)
        except StopIteration:
            print("Number of FASTA sequences not "
                  "equal to number of quality "
                  "sequences", file=sys.stderr)
            return
        if not qual_id == fasta_seq.identifier:
            print("Mismatched FASTA entry and quality score entry:\n"
                  "   FASTA ID: %s\n"
                  "   Quality score ID: %s\n" % (fasta_seq.identifier,qual_id)
                  ,file=sys.stderr)

        elif len(scores) != len(fasta_seq.sequence):
            print("FASTA sequence: ", fasta_seq.identifier,"\n",
                  "is not same length as quality score sequence"
                  " with same identifier. Omitting sequence."
                  ,file=sys.stderr)
        else: yield fastq_seq(fasta_seq.comment,fasta_seq.identifier,fasta_seq.sequence,scores)
#At this point, we've read all the quality score entries, so to make sure the FASTA file had the
#same number of entries, try to read one more entry from the FASTA file, and print an error if
#there are still entries left to get (as this means there were more entries in the FASTA file 
#than were in the quality score file.
    try:
        fasta_seq = next(fasta_seqs)
    except StopIteration:
        pass
    else:
        print("Number of FASTA sequences not "
              "equal to number of quality "
              "sequences", file=sys.stderr)
        return

