#!/bin/bash
set -e
set -x
set -o pipefail

# generate a chromosome renaming file
echo > "/beagle/ref/chrom_rename.txt" <<EOF
chr1 1
chr2 2
chr3 3
chr4 4
chr5 5
chr6 6
chr7 7
chr8 8
chr9 9
chr10 10
chr11 11
chr12 12
chr13 13
chr14 14
chr15 15
chr16 16
chr17 17
chr18 18
chr19 19
chr20 20
chr21 21
chr22 22
chrX X
EOF

# do each chromosome details separately in parallel
# this keeps all chromosomes in one docker image layer, and improves speed
# use xargs to limit concurrent operations
# nproc gets number of processors avaliable, including docker cpu set limits
xargs --max-procs `nproc` --max-args 2 /bin/bash /beagle/ref_41_chrom.sh <<EOF
01
chr1
02
chr2
03
chr3
04
chr4
05
chr5
06
chr6
07
chr7
08
chr8
09
chr9
10
chr10
11
chr11
12
chr12
13
chr13
14
chr14
15
chr15
16
chr16
17
chr17
18
chr18
19
chr19
20
chr20
21
chr21
22
chr22
EOF