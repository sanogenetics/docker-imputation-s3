#!/bin/bash
set -e
set -x
set -o pipefail

# downloads reference 1000Genomes files and pre-processes them as necessary
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=1
# directory to output to
TARGET=$1
mkdir -p "${TARGET}"

# get maps
# get binaries
# run in parallel
# TODO only get one we want, each is only 23-29MB
if [ ! -f "${TARGET}/plink.GRCh38.map.zip" ]; then
  curl -s "http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip" > "${TARGET}/plink.GRCh38.map.zip"
  rm -rf "${TARGET}/map38/"
  mkdir -p "${TARGET}/map38/"
  unzip "${TARGET}/plink.GRCh38.map.zip" -d "${TARGET}/map38/"
fi
if [ ! -f "${TARGET}/plink.GRCh37.map.zip" ]; then
  curl -s "http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip" > "${TARGET}/plink.GRCh37.map.zip"
  rm -rf "${TARGET}/map37/"
  mkdir -p "${TARGET}/map37/"
  unzip "${TARGET}/plink.GRCh37.map.zip" -d "${TARGET}/map37/"
fi
if [ ! -f "${TARGET}/plink.GRCh36.map.zip" ]; then
  curl -s "http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh36.map.zip" > "${TARGET}/plink.GRCh36.map.zip"
  rm -rf "${TARGET}/map36/"
  mkdir -p "${TARGET}/map36/"
  unzip "${TARGET}/plink.GRCh36.map.zip" -d "${TARGET}/map36/"
fi
if [ ! -f "bin/beagle.jar" ]; then
  curl -s "https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar" > "bin/beagle.jar" &
fi
if [ ! -f "bin/conform-gt.jar" ]; then
  curl -s "https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar" > "bin/conform-gt.jar" &
fi
wait

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
#chroms_num=( "1" )
#chroms_z=( "01" )
#chroms_long=( "chr1" )

# generate a chromosome renaming file
for i in "${!chroms_long[@]}"; do
  echo ${chroms_long[i]} ${chroms_num[i]}
done >> "${TARGET}/chrom_rename.txt"

# do each chromosome details separately in parallel
# TODO use xargs to limit concurrent operations
for i in "${!chroms_long[@]}"; do
  /bin/bash ref_2_chr.sh ${TARGET} ${chroms_num[i]} ${chroms_z[i]} ${chroms_long[i]}
done
wait

