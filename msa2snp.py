#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  msa2snp.py
#  
#  Copyright 2019 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	seqnames = []
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line: # skip blank lines
				if line.startswith(">"):
					sequence_name = line.lstrip("> ").split()[0] # left strip > or space, so " > abc edf" will be "abc edf", then split by space to get "abc"
					fasta[sequence_name] = ""
					seqnames.append(sequence_name)
				else:
					fasta[sequence_name] += line.replace(" ", "") # remove spaces in case
	return fasta, seqnames

import sys
aln_file = sys.argv[1] # algnment file, fasta format, first one is the reference

fasta, seqnames = get_fasta(aln_file)
refid = seqnames[0] # reference name
refseq = fasta[refid] # ref seq
refseq_nogap = refseq.replace("-","")

n = -1
seq2 = {} # new dictionary of transformed sequences with no-gap postions as keys

for i in range(len(refseq)): # i is the position of seq poition
	if refseq[i] != "-":
		n += 1
	templist = []
	for j in seqnames: # j is seq name
		seq = fasta[j]
		templist.append(seq[i])
	if n not in seq2:
		seq2[n] = templist
	else:
		seq2[n] = [s1 + s2 for s1, s2 in zip(seq2[n], templist)]
	

# output
outlist = ["Pos\tRef\talt\t" + "\t".join(seqnames)]

for k in range(len(refseq_nogap)):
	seq = [w.replace('-', '') if len(w) > 1 else w for w in seq2[k] ] # remove "-" in alleles in "A--"
	alleles = set(seq) # set: unique values
	if len(alleles) != 1:
		ref_allele = refseq_nogap[k] # string
		alt_allele_set = alleles - set(ref_allele) # set
		alt_allele = "/".join(alt_allele_set) # string
		calls = "\t".join(seq) # SNP calls in each sequence
		outlist.append("\t".join([str(k+1), ref_allele, alt_allele, calls]))

for i in outlist:
	print i		

