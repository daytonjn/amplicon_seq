# Amplicon-seq processing & variant calling

**Description:** Early development scripts for processing multiplexed amplicon-seq, based on barcodes described in Liu et al. (2021). Updated scripts are within period_ms repository. 

**Purpose:** In order to measure the effects of allelic variation & gene knockout on seasonal/daily phenotypes in *O. nubilalis*, an efficient genotyping procedure is necessary. Multiplexed Amplicon-seq data are processed to determined genotypes (VCF) for both gene-editing target loci & known candidate modifer loci.

**Amplicon-seq library preparation workflow:**
<p align="center">
  <img src="https://github.com/user-attachments/assets/b0446396-a421-46d5-8bc3-d1d4e613e3a8" alt="barcode_pcr graphic"/>
</p>

**Workflow:**
  1) Unzip & trim adapters from .fastq output (trim_galore)
  2) Demultiplex .fastq file using .fasta containing barcodes_fwd.fasta & barcodes_rev.fasta (cutadapt)
  3) Rename demultiplexed files, merge paired-end reads (NGmerge), & remove sequences < 110 bp (seqtk)
  4) Align sequence files to reference genome (hisat2) & output BAM alignment (samtools)
  5) Call variants with GATK pipeline (gatk)
  6) Filter variants (hard-filters) & output separate VCFs for SNPs & indels (CRISPR/Cas9-edit)
