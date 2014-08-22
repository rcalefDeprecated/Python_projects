#! /usr/bin/env python3.3
import sys
import re
import collections
import operator
import argparse


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

def get_next_sequence(input):	
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
	return

def main(args):
#	order = parse_arg(args)
#	words = collections.defaultdict(int)
	for sequence in get_next_sequence(sys.stdin):
		sys.stdout.write("\n>%s   %s" % (sequence.identifier,sequence.comment))
		for i in range(len(sequence.sequence)):
			if i % 80 == 0: print
			sys.stdout.write(sequence.sequence[i])
                print
#	print_dict(words, order)


def read_word(input):
	for line in input:
		for word in re.findall(r"[A-Za-z']*", line):
			if word == "": continue
			else: yield word

if __name__ == "__main__":
	sys.exit(main(sys.argv))
