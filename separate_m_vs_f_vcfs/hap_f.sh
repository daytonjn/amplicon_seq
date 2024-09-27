#!/bin/bash
#SBATCH -J ID_VAR #job name
#SBATCH --time=0-1:30:00 
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH --output=genome_index.%j.out
#SBATCH --error=genome_index.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacob.dayton@tufts.edu

module load anaconda/2021.11
source activate ampliconez

module load gatk4

REFERENCE_GENOME=/cluster/tufts/dopmanlab/Jacob/ampliconez/reference_ampliconez3.fasta
BASE_NAMES_FILE="./female_seqs.txt"
BAM_DIR="./mapped"
OUTPUT_DIR="./vcf_female"
mkdir -p $OUTPUT_DIR

# Loop through each base name in the text file
while IFS= read -r BASE_NAME; do
    # Define the input BAM file and output GVCF file paths
    BAM_FILE="${BAM_DIR}/${BASE_NAME}.sorted.bam"
    GVCF_FILE="${OUTPUT_DIR}/${BASE_NAME}.g.vcf"

    # Run GATK HaplotypeCaller to generate GVCF
    gatk HaplotypeCaller \
        -R $REFERENCE_GENOME \
        -I $BAM_FILE \
        -O $GVCF_FILE \
        -ERC GVCF \
        --sample-ploidy 1

    echo "Generated GVCF for ${BASE_NAME}"
done < $BASE_NAMES_FILE