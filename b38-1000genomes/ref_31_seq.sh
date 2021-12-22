#!/bin/bash
set -e
set -x
set -o pipefail

# downloads reference 1000Genomes files and pre-processes them as necessary
# https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw?step=1
CHROM_NUM=$1
CHROM_Z=$2

# download each chromosome reference if not already done
# name to alpha sorting pattern filename
# need the fasta sequence of this chromosome
# download genome reference
# rename it to have chr prefixed to number e.g. chr1
# recompress into block gzip format on the way
curl -s "ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${CHROM_NUM}.fa.gz" | gunzip -c | sed 's/^>\(.*\)$/>chr\1/' | bgzip > "/beagle/ref/${CHROM_Z}.fasta"
# build an index
samtools faidx "/beagle/ref/${CHROM_Z}.fasta"
