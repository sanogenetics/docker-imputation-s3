#!/bin/bash
set -e
set -x
set -o pipefail

# do each chromosome details separately in parallel
# this keeps all chromosomes in one docker image layer, and improves speed
# use xargs to limit concurrent operations to number of cpus
# nproc gets number of processors avaliable, including docker cpu set limits
xargs --max-procs `nproc` --max-args 2 /bin/bash "/beagle/ref_31_seq.sh" <<EOF
1
01
2
02
3
03
4
04
5
05
6
06
7
07
8
08
9
09
10
10
11
11
12
12
13
13
14
14
15
15
16
16
17
17
18
18
19
19
20
20
21
21
22
22
EOF
