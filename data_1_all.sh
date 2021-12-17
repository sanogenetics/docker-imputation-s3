#!/bin/bash
set -e
set -x
set -o pipefail

# prepares a given VCF file for imputation
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=2

VCF=$1
WORKDIR=$2
REF=$3
OUT=$4

# wipe the work area
rm -rf "${WORKDIR}"
mkdir -p "${WORKDIR}"

# various lists of chromosome names
chroms_num=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" )
chroms_z=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" )
chroms_long=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" )
# excluding X
chroms_num=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" )
chroms_z=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" )
chroms_long=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" )

# use only 9 & 22 for testing
# takes ~20m from AWS
# 9 is the shortest 0-padded chromosome
#chroms_num=( "9" "22" )
#chroms_z=( "09" "22" )
#chroms_long=( "chr9" "chr22" )

# split each chromosome
# use xargs to limit concurrent operations otherwise out of memory
rm -f "${WORKDIR}/data_2_chr.sh"
for i in "${!chroms_long[@]}"; do
  echo /bin/bash data_2_chr.sh "${VCF}" ${chroms_num[i]} ${chroms_z[i]} ${chroms_long[i]} "${WORKDIR}" "${REF}" >> "${WORKDIR}/data_2_chr.sh"
done
xargs --max-args 1 --arg-file "${WORKDIR}/data_2_chr.sh" --max-procs 8 --replace bash -c "{}"


# combine output files together
if [ -n `dirname "${OUT}"` ]; then
  # output file has a directory component
  # ensure dir exists
  mkdir -p `dirname "${OUT}"`
fi
bcftools concat --no-version "${WORKDIR}"/reheader.*.vcf.gz -Oz -o "${OUT}"