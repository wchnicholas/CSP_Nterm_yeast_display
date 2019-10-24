#!/usr/bin/python
import os
import sys
import glob
import string
import operator 
from Bio import SeqIO
from string import atof
from itertools import imap
from collections import Counter, defaultdict

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def seq_to_barcode(seq):
  barcode = ''
  for n in range(len(seq)):
    if n%3 == 0:
      barcode+=seq[n]
  return barcode

def refseq_to_barcode(refseq_dict):
  barcode_dict = defaultdict(list)
  for ID in sorted(refseq_dict.keys()):
    refseq   = refseq_dict[ID]
    barcode  = seq_to_barcode(refseq)
    barcode1 = barcode.replace('K','G')
    barcode2 = barcode.replace('K','T')
    barcode_dict[barcode1].append(ID)
    barcode_dict[barcode2].append(ID)
  return barcode_dict

def callsub(mutpep,refseq):
  shift = 0
  haplo = []
  assert(len(mutpep)==len(refseq))
  for n in range(len(mutpep)):
    pos = n+shift
    if refseq[n]!=mutpep[n]:
       haplo.append(refseq[n]+str(pos)+mutpep[n])
  if len(haplo) == 0: return 'WT'
  else: return '-'.join(haplo)

def ReadingRefSeq(refseq_file):
  refseq_dict = {}
  records = SeqIO.parse(refseq_file,"fasta")
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    refseq_dict[ID] = seq
  return refseq_dict

def Read2Mut2Count(R1_file, Count_dict, Sample, refseq_dict, barcode_dict):
  print "Reading %s" % R1_file
  R2_file = R1_file.replace('_R1','_R2')
  FPrimerlength = 25
  RPrimerlength = 25
  roilength     = 218
  R1records = SeqIO.parse(R1_file,"fastq")
  R2records = SeqIO.parse(R2_file,"fastq")
  variants = []
  WTpep    = translation(refseq_dict['0'][1:-1])
  for R1record in R1records:
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi = R1seq[FPrimerlength:FPrimerlength+roilength]
    R2roi = R2seq[RPrimerlength:RPrimerlength+roilength]
    if 'N' in R1roi or 'N' in R2roi: continue
    if R1roi == rc(R2roi):
      barcode = seq_to_barcode(R1roi)
      R1pep   = translation(R1roi[1:-1])
      if barcode not in barcode_dict.keys(): continue
      mut_positions = map(int,barcode_dict[barcode])
      mut_haplo = []
      for mut_position in mut_positions:
        WTaa  = WTpep[mut_position-1]
        mutaa = R1pep[mut_position-1]
        if WTaa != mutaa:
          mut = WTaa+str(mut_position)+mutaa
          mut_haplo.append(mut)
      if len(mut_haplo) == 0:
        sub = 'WT'
      elif len(mut_haplo) == 1:
        sub = mut_haplo[0]
      elif len(mut_haplo) > 1:
        continue
      Count_dict[Sample][sub] += 1
  return Count_dict

def Output(Count_dict, Sample_outfile, position_offset):
  Samples = Count_dict.keys()
  outfile = open(Sample_outfile,'w')
  print "Compiling results into files with prefix: %s" % Sample_outfile
  Muts    = list(set([mut for S in Count_dict.keys() for mut in Count_dict[S].keys()]))
  outfile.write("\t".join(map(str,['Mut','Sample','Count']))+"\n")
  for Sample in Samples:
    for Mut in Muts:
      count = Count_dict[Sample][Mut]
      Mut = Mut[0]+str(int(Mut[1:-1])+position_offset)+Mut[-1] if Mut != 'WT' else 'WT'
      outfile.write("\t".join(map(str,[Mut, Sample, count]))+"\n")
  outfile.close()

def main():
  position_offset = 25
  refseq_file    = 'Fasta/ref_seqs.fa'
  Sample_outfile = 'result/nterm_CSP_sub_count.tsv'
  refseq_dict    = ReadingRefSeq(refseq_file)
  barcode_dict   = refseq_to_barcode(refseq_dict)
  R1_files = glob.glob('fastq/*_R1*.fastq')
  Count_dict  = {}
  for R1_file in sorted(R1_files):
    Sample = R1_file.rsplit('/')[1].rsplit('_')[0]
    if Sample not in Count_dict.keys(): Count_dict[Sample] = defaultdict(int)
    Count_dict = Read2Mut2Count(R1_file, Count_dict, Sample, refseq_dict, barcode_dict)
  Output(Count_dict, Sample_outfile, position_offset)

if __name__ == "__main__":
  main()
