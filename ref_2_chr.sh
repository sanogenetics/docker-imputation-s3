#!/bin/bash
set -e
set -x
set -o pipefail

# downloads reference 1000Genomes files and pre-processes them as necessary
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=1
# directory to output to
TARGET=$1
CHROM_NUM=$2
CHROM_Z=$3
CHROM_LONG=$4

# download each chromosome reference if not already done
# name to alpha sorting pattern filename
# need the fasta sequence of this chromosome
if [ ! -f "${TARGET}/${CHROM_Z}.fasta" ]; then
  # download genome reference
  # recompress into block gzip format on the way
  curl -s "ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${CHROM_NUM}.fa.gz" | gunzip -c | bgzip > "${TARGET}/${CHROM_Z}.fasta"
fi

if [ ! -f "${TARGET}/${CHROM_Z}.fasta.fai" ]; then
  # build an index
  samtools faidx "${TARGET}/${CHROM_Z}.fasta"
fi

if [ ! -f "${TARGET}/${CHROM_Z}.vcf.gz" ]; then
  # download it
  # rename chromosome to only numeric (no chr prefix)
  # remove variants seen in only 1 or 2 samples
  # split multiallelic sites to biallelic records
  # keep only SNPs and INDELs
  # align to reference so REF and ALT in the shortest form and that REF matches reference, additionally remove duplicate variants
  # after alignment, remove multiallelic records since these are formed during the alignment if the REF does not match the reference
  # remove sites containing missing data
  curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.${CHROM_LONG}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz" | \
    bcftools annotate --no-version --rename-chrs "${TARGET}/chrom_rename.txt" -Ou | \
    bcftools view --no-version -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou | \
    bcftools norm --no-version -m -any -Ou | \
    bcftools view --no-version -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm --no-version -f "${TARGET}/${CHROM_Z}.fasta" -d none -Ou | \
    bcftools view --no-version -m 2 -M 2 -Ou | \
    bcftools view --no-version -g ^miss -Oz -o "${TARGET}/${CHROM_Z}.vcf.gz"
fi

# TODO handle X chromosome ploidy to always diploid

# TODO remove RSID duplicates (if any)

# TODO convert to bref (optional)
# https://faculty.washington.edu/browning/beagle/bref3.28Jun21.220.jar


