Beagle imputation
=================

Takes input and output in VCF format.

input must:

- have chromosomes named either with or without chr prefix (e.g. 1, chr1)

output will:

- have chromosomes named with chr prefix (e.g. chr1)

Based on process described in https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw 

TODO
----
X chromosome
b37 
b36
population specifc reference panels
more reference panel sources

scripts
-------

`buildall.sh`  
builds all docker images and tags them as `b38-1000genomes-latest` etc

Docker recap
------------

`docker system df -v`
List current images and sizes

`docker system prune`
Cleanup danglers

`docker system prune -a`
Cleanup everything!