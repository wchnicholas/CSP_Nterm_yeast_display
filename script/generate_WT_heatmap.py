#!/usr/bin/python
import os
import sys

def aa_index_mapper(aa):
  aa_list = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  aa_list.reverse()
  return aa_list.index(aa)+1

def generate_aa_index(seq, starting_pos, outfile):
  print "writing: %s" % outfile
  count = 0
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['resi','x','y'])+"\n")
  for aa in seq:
    count += 1
    aa_index = aa_index_mapper(aa)
    outfile.write("\t".join(map(str,[aa+str(starting_pos), count, aa_index]))+"\n")
    starting_pos += 1
  outfile.close()

def main():
  fasta_file = "Fasta/WT_pep.fa"
  outfile    = "data/WT_heatmap.tsv"
  starting_pos = 26
  seq = open(fasta_file,'r').readlines()[1].rstrip()
  generate_aa_index(seq, starting_pos, outfile)

if __name__ == "__main__":
  main()
