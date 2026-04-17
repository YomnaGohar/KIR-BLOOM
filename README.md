![KIR*bloom logo](logo2.jpeg)
This repository contains a newly developed method for accurate KIR allele genotyping from short-read sequencing data

The pipeline accepts fasta, bam, cram file.
for bam or cram you need to specify the reference it was mapped to or you have to provide a bed file
of the kir regions.
for each type of files you need to use a config file.
make sure that for the fastq files that they are paired and pair 1 ends with /1 and pair2 ends with /2 otherwise bwa mem will give you an error.