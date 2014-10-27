#! /usr/bin/env python
import string
from collections import Counter

def count_kmers(sequence,order,start,stop,alphabet):
	alphabet = alphabet + start + stop
	alphabet = set(alphabet)
	counts = Counter()
	k = order + 1
	if k < 3: 
		seq = string.upper(start + sequence + stop)
	else: 
		seq = string.upper(start*order + sequence + stop*order)
	for kmer in kmers_from_sequence(seq,alphabet,k):
		counts[kmer] += 1
	return counts


def simple_kmers_from_sequence(sequence,k):
	for start in range(len(sequence) - k + 1):
		yield sequence[start:start+k]

def kmers_from_sequence(sequence,alphabet,k):
	upper_bound = len(sequence)-k+1
	for itor in range(upper_bound):
		if sequence[itor] not in alphabet: continue
		kmer=""
		char_check=itor
#Want to build up next kmer, skipping characters not in alphabet
		while len(kmer) < k:
#If char_check runs off end of sequence, then no more kmeers, return
			if char_check == len(sequence): return
			if sequence[char_check] in alphabet: kmer += sequence[char_check]
			char_check += 1
		yield kmer

def get_permutations(alphabet,length):
	alpha_len = len(alphabet)
	indices=[0 for x in range(length)]
	element=""
	more = True
	alphabet=list(alphabet)
	while(more):
		for index in indices:
			element += alphabet[index]
		yield element
		element=""
		for i in range(length):
			indices[i] += 1
			if indices[i] != alpha_len:
				break
			else:
				if i == length-1:
					more=False
				indices[i]=0

def get_kmers(start,stop,alphabet,k):
	alpha_len = len(alphabet)
#First get kmers with start characters
	if k == 1:
		yield start
	else:
		for i in range(1,k):
			for mer in get_permutations(alphabet,k-i):
				yield (i*start) + mer
#Then get kmers without start or stop characters
	for kmer in get_permutations(alphabet,k):
		yield kmer
#Finally get kmers with end characters
	if k == 1:
		yield stop
	else:
		for i in range(1,k):
			for mer in get_permutations(alphabet,k-i):
				yield mer + (i*stop)
