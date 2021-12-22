#!/bin/bash
set -e
set -x
set -o pipefail

mkdir -p "/beagle/bin"

# get binaries
curl -s "https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar" > "/beagle/bin/beagle.jar" 
curl -s "https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar" > "/beagle/bin/conform-gt.jar"
curl -s "https://faculty.washington.edu/browning/beagle/bref3.28Jun21.220.jar" > "/beagle/bin/bref3.jar"

