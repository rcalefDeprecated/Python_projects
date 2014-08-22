#! /usr/bin/env python3.3
from __future__ import print_function
import sys
import re
import collections
import operator
import argparse

allowed_bases=['A','C','G','T','N','U','K','S','Y','M','W','R','B','D','H','V','-']

class fastq_seq:
	def __init__(self,comment,identifier,sequence,qual_scores):
		self.comment = comment
		self.identifier = identifier
		self.sequence = sequence
		self.scores=qual_scores

class fasta_seq:
	def __init__(self,comment,identifier,sequence):
		self.comment = comment
		self.identifier = identifier
		self.sequence = sequence

def parse_arg(args):
	argparser = argparse.ArgumentParser(description = ("Reads input from std"
		"in, splits input into words containing alphabetic char"
		"acters, and returns pairs of words along with their"
		" number of occurences in the input"))
	asc_or_desc = argparser.add_mutually_exclusive_group(
		required=True)
	asc_or_desc.add_argument('--ascend', action = 'append_const',
		const = 2, dest = 'order', help = ("Returns output "
		"sorted in ascending order by counts of occurrences."))
	asc_or_desc.add_argument('--descend', action = 'append_const',
		const = 4, dest = 'order', help = ("Returns output "
		"sorted in descending order by counts of occurrences."))
	argparser.add_argument('--alpha', action = 'append_const',
		const = 1, dest = 'order', help =("Returns output "
		"sorted alphabetically by word. If specified along"
		"with --descend, output is sorted in reverse alpha"
		"betic order. Specifying --alpha and --ascend gives"
		"the same results as only --alpha."))
	options = argparser.parse_args()
	print(vars(options))
	order = 0
	for num in options.order:
		order += num
	print(order)
	return order
	
def print_dict(words, order):
	first_sort = 1
	second_sort = 0
	rvrse = False
	if 0b1 & order:
		first_sort = 0
		second_sort = 1
	if 0b10 & order:
		pass
	if 0b100 & order:
		rvrse = True
	sorted_words = sorted(words.items(),
		key=operator.itemgetter(first_sort,second_sort), reverse=rvrse)
	for entry in sorted_words:
		print(entry[0], entry[1])

def read_fastq(input):	
	bases_read=scores_read=0
	in_scores=in_sequence=False
	match=seq_id=comment=sequence=""
	scores=[]
#Now iterate over input until a new sequence is found
	for line in input:
#If first character of line is '>', then we've encountered the 
# start of a new sequence entry, need to set this line to be
#curr_line and return the current sequence.
                if in_scores and bases_read==scores_read:
			yield fastq_seq(comment,seq_id,sequence,scores)
			in_sequence=False
			in_scores=False
                if not in_sequence:
                        if line[0] !='@':
                                continue
                        match = re.search("[\s,]",line)
                        seq_id=line[1:match.start(0)]
                        comment=line[match.start(0):]
                        sequence=""
                        scores=[]
                        in_sequence = True
                        bases_read=scores_read=0
                        continue

		if line.isspace(): continue
		if not in_scores and line[0] == '+':
			in_scores = True
			continue
#Check if valid character
		for char in line:
#If white space, continue
			if char.isspace(): continue
			if not in_scores:
				if not char.upper()in allowed_bases: continue
				else:
					sequence += char
					bases_read += 1
			else:
				score=ord(char)-33
				if score < 0: continue
				scores.append(score)
				scores_read +=1
#If we reach the end of file, then return last entry
	yield fastq_seq(comment,seq_id,sequence,scores)

def read_fasta(input):
#If we haven't found a sequence entry yet, find the beginning of the first one
        for line in input:
#Read through any whitespace or comments until beginning of sequence
                if line[0] != '>': continue
                else: break
#Split first line into sequence ID and comment.
        match = re.search("[\s,]",line)
        seq_id=line[1:match.start(0)]
        comment=line[match.start(0):]
        sequence=""
#Now iterate over input until a new sequence is found
        for line in input:
#If first character of line is '>', then we've encountered the
# start of a new sequence entry, need to set this line to be
#curr_line and return the current sequence.
                if line[0] == '>':
                        yield fasta_seq(comment,seq_id,sequence)
                        match = re.search("[\s,]",line)
                        seq_id=line[1:match.start(0)]
                        comment=line[match.start(0):]
                        sequence=""
                        continue
                if line.isspace(): continue
#Check if valid character
                for char in line:
#If white space, continue
                        if char.isspace(): continue
                        if char.isalpha():
#Only disallowed alphabet letters in a FASTA sequence are J or O
                                if char.upper() == 'J' or char.upper() == 'O': continue
                                else: sequence += char
#If not a letter, only allowed symbols are * or -
                        elif char == '*' or char == '-': sequence += char
#If we reach the end of file, then return last entry
        yield fasta_seq(comment,seq_id,sequence)

def read_fasta_with_quality(fasta_input, qual_input):
#First get FASTA reader and initialize variables used throughout loop
	fasta_seqs=read_fasta(fasta_input)
	in_score=False
	scores=[]
	for line in qual_input:
#If currently reading a score, and find a line beginning with '>'
#then get the next FASTA sequence, and yield along with scores
#as a new FASTQ sequence
		if in_score and line[0]=='>':
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
				yield None
#Else yield the new sequence and reset in_score to indicate we need
#to start a new sequence.
			else: yield fastq_seq(fasta_seq.comment,fasta_seq.identifier,fasta_seq.sequence,scores)
			in_score=False
                if not in_score:
                        if line[0] != '>': continue
                        match = re.search("[\s,]",line)
                        seq_id=line[1:match.start(0)]
                        scores=[]
                        in_score=True
                        continue
#If we're reading a block (in_score==True) and the line does not begin
# with '>', then split the line into quality scores, check each to
#make sure they're a valid digit, and append to current score sequence
		for num in line.split():
			if(num.isdigit()):
				score=int(num)
				if score >= 0: scores.append(score)
#When we reach end of file, we still have to yield the final sequence,
#hence the duplication of code below
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
		yield None
	else: yield fastq_seq(fasta_seq.comment,fasta_seq.identifier,fasta_seq.sequence,scores)
	


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
#	print_dict(words, order)


def read_word(input):
	for line in input:
		for word in re.findall(r"[A-Za-z']*", line):
			if word == "": continue
			else: yield word

if __name__ == "__main__":
	sys.exit(main(sys.argv))
