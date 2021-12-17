#!/bin/bash
set -e
set -x
set -o pipefail

# prepares a given VCF file for imputation
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=2

VCF=$1
CHROM_NUM=$2
CHROM_Z=$3
CHROM_LONG=$4
WORKDIR=$5
REF=$6

# keep only wanted chromosomes
# align the alleles to the reference genome, 
# keep only biallelic records
bcftools view --no-version -t ${CHROM_NUM} $VCF -Ou | \
  bcftools norm --no-version -f "${REF}/${CHROM_Z}.fasta" -c ws -Ou | \
  bcftools view --no-version -m 2 -M 2 -Ou | \
  bcftools norm --no-version -d none -Oz -o "${WORKDIR}/clean.${CHROM_Z}.vcf.gz"

# check with conform-gt that everything matches
# match=POS use position to pair reference and sample, not ID
rm -f "${WORKDIR}/conform-gt.${CHROM_Z}.vcf.gz"   # remove any old file
/usr/bin/time java -jar bin/conform-gt.jar \
  match=POS \
  chrom=${CHROM_NUM} \
  gt="${WORKDIR}/clean.${CHROM_Z}.vcf.gz" \
  ref="${REF}/${CHROM_Z}.vcf.gz" \
  out="${WORKDIR}/conform-gt.${CHROM_Z}"

# run the imputation
rm -f "${WORKDIR}/beagle.${CHROM_Z}.vcf.gz"   # remove any old file
/usr/bin/time java -jar bin/beagle.jar \
  gt="${WORKDIR}/conform-gt.${CHROM_Z}.vcf.gz" \
  ref="${REF}/${CHROM_Z}.vcf.gz" \
  map="${REF}/map38/plink.${CHROM_LONG}.GRCh38.map" \
  nthreads=16 \
  impute=true \
  out="${WORKDIR}/beagle.${CHROM_Z}"

# after imputation, beagle leaves out an necessary header line
# ##contig=<ID=9>
rm -f "${WORKDIR}/header.${CHROM_Z}.txt"
bcftools view -h "${WORKDIR}/beagle.${CHROM_Z}.vcf.gz" | head -n -1 > "${WORKDIR}/header.${CHROM_Z}.txt"
echo "##contig=<ID=${CHROM_NUM}>" >> "${WORKDIR}/header.${CHROM_Z}.txt"
bcftools view -h "${WORKDIR}/beagle.${CHROM_Z}.vcf.gz" | tail -n 1 >> "${WORKDIR}/header.${CHROM_Z}.txt"
bcftools reheader -h "${WORKDIR}/header.${CHROM_Z}.txt" -o "${WORKDIR}/reheader.${CHROM_Z}.vcf.gz" "${WORKDIR}/beagle.${CHROM_Z}.vcf.gz"