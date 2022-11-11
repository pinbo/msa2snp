#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  msa2vcf.py
#  
#  Copyright 2022 Junli Zhang <zhjl86@gmail.com>
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
# example: python msa2vcf.py msa-example.fa > snp.vcf

usage = '''
msa2vcf.py <alinged-sequences.fa>

alinged-sequences.fa: an aligned fasta file:
        1) 1st sequence is the reference;
        2) flat ends (equal length for all sequences)
        3) NO "-" in the beginning of the alignment for any sequences

Output: the vcf output will be printed in stdout, and a SNP table will be write to disk named "snptable.txt"

Example: python msa2vcf.py msa-example.fa > snp.vcf
'''
version = "0.1.0"

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

# vcf header
vcfheader = '''##fileformat=VCFv4.2
##contig=<ID=mychrom,length=chromlength>
##fileDate=todaydate
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##source=msa2snp.py
##msa2vcf.pyVersion="myversion"
##msa2vcf.pyCmd="mycmd"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'''
vcfheader = vcfheader.replace("myversion", version)

# add date
from datetime import date
today = date.today()
d1 = today.strftime("%Y%m%d")
vcfheader = vcfheader.replace("todaydate", d1)

# get arguments
import sys

if len(sys.argv) != 2:
	print(usage)
	sys.exit(0)

aln_file = sys.argv[1] # algnment file, fasta format, first one is the reference

fasta, seqnames = get_fasta(aln_file)
refid = seqnames[0] # reference name
refseq = fasta[refid] # ref seq
refseq_nogap = refseq.replace("-","")
refseq_lengh = len(refseq_nogap) # reference length

# update header
vcfheader = vcfheader.replace("mychrom", refid)
vcfheader = vcfheader.replace("chromlength", str(refseq_lengh))
vcfheader = vcfheader.replace("mycmd", " ".join(sys.argv))
vcfheader = vcfheader + "\t" + "\t".join(seqnames)
print(vcfheader)
# start processing
n_nogap = -1 # no gap position for reference
n = 0 # no gap position for mutations
seq2 = {} # new dictionary of rotated sequences with no-gap postions as keys

for i in range(len(refseq)): # i is the position of seq poition
	if refseq[i] != "-":
		n_nogap += 1
	templist = []
	nblank = 0
	for j in seqnames: # j is seq name
		seq = fasta[j]
		templist.append(seq[i])
		if seq[i] == "-":
			nblank += 1
	if nblank == 0:
		n = n_nogap
	if n not in seq2:
		seq2[n] = templist
	else:
		seq2[n] = [s1 + s2 for s1, s2 in zip(seq2[n], templist)]

# print(seq2)
# output
outlist = []

for k in range(len(refseq_nogap)):
	if (k in seq2):
		seq = [w.replace('-', '') if len(w) > 1 else w for w in seq2[k] ] # remove "-" in alleles in "A--"
		alleles = set(seq) # set: unique values
		if len(alleles) != 1: # if not monomorphic in this position
			ref_allele = seq[0] #refseq_nogap[k] # string
			alt_allele_set = alleles - set([ref_allele]) # set
			alt_allele = ",".join(alt_allele_set) # string
			outlist.append([str(k+1), ref_allele, alt_allele] + seq) # K+1: starting from 1

# print(outlist)

## format to vcf format
# each line: Pos	Ref	alt	seq1	seq2	seq3
# to vcf line: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	seq1	seq2	seq3
# write vcf to stdout, write snptable to a file
f = open("snptable.txt", "w")
f.write("\t".join(["Pos", "Ref", "alt"] + seqnames) + "\n")
for i in outlist:
	alist = i[1:2] + i[2].split(",") # allele list
	indexlist = [str(alist.index(x))+"/" + str(alist.index(x)) for x in i[3:]]
	j = [refid] + i[0:1] + ["."] + i[1:3] + [".\t.\t.\tGT"] + indexlist
	outline = "\t".join(j)
	print (outline)
	f.write("\t".join(i) + "\n")

f.close()