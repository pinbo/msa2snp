# msa2snp
SNP and indel calls from multiple sequence alignment (MSA)

I need to call SNPs in a multiple sequence alignment from several sanger sequencing sequences, but did not find a simple software that can do that, so I just wrote one.

After you get the SNPs and their locations, then you can use another script [mysnpeffv2.py](https://github.com/pinbo/gene_manual_annotation) to check SNP effects.

# Input

The only input is the MSA file in fasta format. You can prepare the alignment in many software, such as MEGA and Aliview. Then just save them in fasta format. Some format requirement:

1. The first sequence is the reference.

1. Make sure the alignment has flat ends, that is, equal length for all sequences.

1. Make sure to no "-" in the beginning of the alignment, at least for the reference.

# Usage

```sh
./msa2snp.py msa-example.fa > SNPs.txt
```
