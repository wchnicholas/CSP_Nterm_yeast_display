This README describes the scripts used for analyzing the deep mutational scanning of circumsporozoite protein (CSP) N-terminal domain from Plasmodium falciparum using yeast display. Antibody 5D5 was used for selection. Experiment was performed by David Oyen.

### REQUIREMENTS
* [Python](https://www.python.org/) version 2.7
* [R](https://www.r-project.org) version 3.6.1

### INPUT FILES
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA578947](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA578947), should be placed in fastq/ folder. The filenames for read 1 should match those described in [./data/SampleInfo.tsv](./data/SampleInfo.tsv). The filenames for read 2 should be the same as read 1 except "\_R1" is replaced by "\_R2".
* [./data/SampleInfo.tsv](./data/SampleInfo.tsv): Describe the FASTQ filenames of the raw sequencing read
* [./Fasta/WT\_pep.fa](./Fasta/WT\_pep.fa): Protein sequence of the wild type
* [./Fasta/ref\_seqs.fa](./Fasta/ref\_seqs.fa): Reference sequences for mapping

### ANALYSIS PIPELINE
1. [./script/fastq2mut.py](./script/fastq2mut.py): Converts raw reads to variant counts 
    - Input file:
      - Raw sequencing reads in fastq/ folder
    - Output file:
      - [./result/nterm\_CSP\_sub_count.tsv](./result/nterm\_CSP\_sub_count.tsv)
2. [./script/generate_WT_heatmap.py](./script/generate_WT_heatmap.py): Generate a file that describes the location of the wild type residues when plotting heatmap
    - Input file:
      - [./Fasta/WT\_pep.fa](./Fasta/WT\_pep.fa)
    - Output file:
      - [./data/WT_heatmap.tsv](./data/WT_heatmap.tsv)
3. [./script/plot_enrich_heatmap.R](./script/plot_enrich_heatmap.R): Plot the enrichment heatmap that represents the 5D5 selection
    - Input file:
      - [./result/nterm\_CSP\_sub_count.tsv](./result/nterm\_CSP\_sub_count.tsv)
      - [./data/WT_heatmap.tsv](./data/WT_heatmap.tsv)
    - Output file:
      - [./graph/norm_affinity_heatmap.png](./graph/norm_affinity_heatmap.png)
