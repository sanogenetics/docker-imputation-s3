#!/bin/bash
set -e
set -x
set -o pipefail

mkdir -p "/beagle/ref"

# download maps
curl -s "http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip" > "/beagle/ref/plink.GRCh38.map.zip"
# extract individual maps
unzip "/beagle/ref/plink.GRCh38.map.zip" -d "/beagle/ref/"
# cleanup the downloaded file
rm "/beagle/ref/plink.GRCh38.map.zip"

# each map uses numeric chromosome ids (e.g. 8) but we want with chr prefix (e.g. chr8)
for FILE in /beagle/ref/plink.*.GRCh38.map; do
  sed 's/^\(.*\)$/chr\1/' "${FILE}" > "${FILE}.chrom"
  rm "${FILE}"
  mv "${FILE}.chrom" "${FILE}"
done