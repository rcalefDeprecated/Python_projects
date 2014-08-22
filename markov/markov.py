#! /usr/bin/env python
from collections import Counter
import string

def count_kmers(sequence,order,start,stop,alphabet):
	if alphabet is None: alphabet = string.ascii_uppercase
	alphabet = alphabet + start + stop
	alphabet = set(alphabet)
	counts = Counter()
	k = order + 1
	if k < 3: 
		seq = string.upper(start + sequence + stop)
	else: 
		seq = string.upper(start*order + sequence + stop*order)
	upper_bound = len(seq)-order
	for itor in range(upper_bound):
		if seq[itor] not in alphabet: continue
		mer = ""
		char_check=itor
		while len(mer) < k:
			if char_check == len(seq): return counts
			if seq[char_check] in alphabet: mer += seq[char_check]
			char_check += 1
		counts[mer] += 1
	return counts


