#!/bin/bash
set -e
set -x
set -o pipefail

# prepares a given VCF file for imputation
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=2

VCF=$1
CHROM_Z=$2
CHROM_LONG=$3

# keep only wanted chromosomes
# align the alleles to the reference genome
# keep only biallelic records
bcftools annotate --no-version --rename-chrs "/beagle/ref/chrom_rename.txt" "${VCF}" -Ou | \
  bcftools view --no-version -t ${CHROM_LONG} -Ou | \
  bcftools norm --no-version -f "/beagle/ref/${CHROM_Z}.fasta" -c ws -Ou | \
  bcftools view --no-version -m 2 -M 2 -Ou | \
  bcftools norm --no-version -d none -Oz -o "/beagle/wrk/clean.${CHROM_Z}.vcf.gz"

# check with conform-gt that everything matches
# match=POS use position to pair reference and sample, not ID
# wants reference in vcf format
rm -f "/beagle/wrk/conform-gt.${CHROM_Z}.vcf.gz"   # remove any old file
java -jar bin/conform-gt.jar \
  match=POS \
  chrom=${CHROM_LONG} \
  gt="/beagle/wrk/clean.${CHROM_Z}.vcf.gz" \
  ref="/beagle/ref/${CHROM_Z}.vcf.gz" \
  out="/beagle/wrk/conform-gt.${CHROM_Z}"

# run the imputation
# use the VCF that conform-gt produces that matches the reference
# wants reference in bref3 format, can accept vcf but 33% slower
rm -f "/beagle/wrk/beagle.${CHROM_Z}.vcf.gz"   # remove any old file
java -jar bin/beagle.jar \
  gt="/beagle/wrk/conform-gt.${CHROM_Z}.vcf.gz" \
  ref="/beagle/ref/${CHROM_Z}.bref3" \
  map="/beagle/ref/plink.${CHROM_LONG}.GRCh38.map" \
  nthreads=16 \
  impute=true \
  out="/beagle/wrk/beagle.${CHROM_Z}"

# after imputation, beagle leaves out an necessary header line
# ##contig=<ID=9>
rm -f "/beagle/wrk/header.${CHROM_Z}.txt"
bcftools view -h "/beagle/wrk/beagle.${CHROM_Z}.vcf.gz" | head -n -1 > "/beagle/wrk/header.${CHROM_Z}.txt"
echo "##contig=<ID=${CHROM_LONG}>" >> "/beagle/wrk/header.${CHROM_Z}.txt"
bcftools view -h "/beagle/wrk/beagle.${CHROM_Z}.vcf.gz" | tail -n 1 >> "/beagle/wrk/header.${CHROM_Z}.txt"
bcftools reheader -h "/beagle/wrk/header.${CHROM_Z}.txt" -o "/beagle/wrk/reheader.${CHROM_Z}.vcf.gz" "/beagle/wrk/beagle.${CHROM_Z}.vcf.gz"