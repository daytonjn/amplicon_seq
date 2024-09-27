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

#module load anaconda/2021.11
#source activate ampliconez

#module load gatk4

LIBRARY_NAME=2024_dd
REFERENCE_GENOME=/cluster/tufts/dopmanlab/Jacob/ampliconez/reference_ampliconez3.fasta

mkdir vcf_fm
cp ./vcf_male/*vcf* ./vcf_fm/
cp ./vcf_female/*vcf* ./vcf_fm/
cp ./vcf_unknown/*vcf* ./vcf_fm/

LIBRARY_NAME=2024

# Step 3: Create a sample list file
ls ./vcf_fm/*.g.vcf > ./vcf_fm/sample_gvcfs.list

# Step 4: Combine GVCFs into a single VCF
gatk CombineGVCFs \
    -R $REFERENCE_GENOME \
    --variant ./vcf_fm/sample_gvcfs.list \
    -O ./vcf_fm/${LIBRARY_NAME}_combined.g.vcf

# Step 5: Perform joint genotyping
gatk GenotypeGVCFs \
    -R $REFERENCE_GENOME \
    -V ./vcf_fm/${LIBRARY_NAME}_combined.g.vcf \
    -O ./vcf_fm/${LIBRARY_NAME}_combined_variants.vcf


# Step 6: Optional: Filter variants
gatk VariantFiltration \
    -R $REFERENCE_GENOME \
    -V ./vcf_fm/${LIBRARY_NAME}_combined_variants.vcf \
    -O ./vcf_fm/${LIBRARY_NAME}_filtered_variants.vcf \
    --filter-name "my_filter" \
    --filter-expression "QD < 2.0 || QUAL < 30 || MQ < 40.0"
    
grep -v "^#" ./vcf_fm/${LIBRARY_NAME}_filtered_variants.vcf | cut -f 7 | sort | uniq -c

# Step 7: Select only passing variants:

gatk SelectVariants \
   -R $REFERENCE_GENOME \
   -V ./vcf_fm/${LIBRARY_NAME}_filtered_variants.vcf \
   -O ./vcf_fm/${LIBRARY_NAME}_final_indels.vcf \
   --exclude-filtered \
   --select-type INDEL
   
gatk SelectVariants \
   -R $REFERENCE_GENOME \
   -V ./vcf_fm/${LIBRARY_NAME}_filtered_variants.vcf \
   -O ./vcf_fm/${LIBRARY_NAME}_SNPs.vcf \
   --exclude-filtered \
   --select-type SNP

gatk VariantFiltration \
    -V ./vcf_fm/${LIBRARY_NAME}_SNPs.vcf \
    --filter-expression "AF <= 0.15 || AF >= 0.85" \
    --filter-name "uninformativeAF" \
    -O ./vcf_fm/${LIBRARY_NAME}_filtered_SNPs.vcf
    
    gatk SelectVariants \
       -R $REFERENCE_GENOME \
       -V ./vcf_fm/${LIBRARY_NAME}_filtered_SNPs.vcf \
       -O ./vcf_fm/${LIBRARY_NAME}_final_SNPs.vcf \
       --exclude-filtered
       
       grep -v "^#" ./vcf_fm/${LIBRARY_NAME}_final_SNPs.vcf | cut -f 7 | sort | uniq -c
