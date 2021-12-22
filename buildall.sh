#!/bin/bash
set -e
set -x
set -o pipefail

/usr/bin/time docker build --tag sano-impute:b38-1000genomes-latest ./b38-1000genomes 
