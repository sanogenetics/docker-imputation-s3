#!/bin/bash
set -e
set -x
set -o pipefail

INPUT=$1
OUTDIR=$2
N=5

# print a common header
mkdir -p "${OUTDIR}"
bcftools view --header-only "${INPUT}" > "${OUTDIR}/header.txt"
gzip -f "${OUTDIR}/header.txt"
echo "wrote header file"

# separate the header 
# shuffle the lines
#   gshuf is gnu shuf on mac via brew
# split into separate files
#   gsplit is gnu split on mac via brew
mkdir -p "${OUTDIR}/parts"
bcftools view --no-header "${INPUT}" | shuf | split -n "r/${N}" --filter='gzip > $FILE.gz' --numeric-suffixes=1 - "${OUTDIR}/parts/chunk_"

for i in `seq -f "%02g" ${N}`
do 
  TOCAT="${OUTDIR}/header.txt.gz"
  for j in `seq -f "%02g" ${N}`
  do
    if [[ $i -ne $j ]]
    then
        TOCAT="$TOCAT ${OUTDIR}/parts/chunk_${j}.gz"
    fi
  done
  echo $TOCAT

  
  gunzip -c ${TOCAT} | bgzip > "${OUTDIR}/p${i}.unsorted.vcf.gz"
  bcftools sort "${OUTDIR}/p${i}.unsorted.vcf.gz" -Oz -o "${OUTDIR}/p${i}.vcf.gz"
  rm "${OUTDIR}/p${i}.unsorted.vcf.gz"
done

# cleanup
rm -r "${OUTDIR}/parts"
rm "${OUTDIR}/header.txt.gz"
