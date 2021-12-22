#!/bin/bash
set -e
set -x
set -o pipefail

# prerequisite package instalation
# assumes debian base
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get -y upgrade
apt-get -y install --no-install-recommends curl unzip gzip bcftools tabix samtools
apt-get clean