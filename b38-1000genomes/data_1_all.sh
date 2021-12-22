#!/bin/bash
set -e
set -x
set -o pipefail

# prepares a given VCF file for imputation
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=2

VCF=$1
OUT=$2

# wipe the work area
# this could be in the image, or could be a bind mount
mkdir -p /beagle/wrk
rm -rf /beagle/wrk/*

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

# generate a chromosome renaming file
cat  <<EOF > "/beagle/ref/chrom_rename.txt"
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
EOF

# split each chromosome
# use xargs to limit concurrent operations otherwise out of memory
rm -f /beagle/wrk/data_2_chr.sh
for i in "${!chroms_long[@]}"; do
  echo /bin/bash /beagle/data_2_chr.sh "${VCF}" ${chroms_z[i]} ${chroms_long[i]} >> /beagle/wrk/data_2_chr.sh
done
xargs --max-args 1 --arg-file /beagle/wrk/data_2_chr.sh --max-procs `nproc` --replace bash -c "{}"
rm -f /beagle/wrk/data_2_chr.sh


# output file has a directory component
if [ -n `dirname "${OUT}"` ]; then
  # ensure dir exists
  mkdir -p `dirname "${OUT}"`
fi
# combine output files together
bcftools concat --no-version /beagle/wrk/reheader.*.vcf.gz -Oz -o "${OUT}"