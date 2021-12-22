Beagle imputation
=================

Takes input and output in VCF format.

input must:

- have chromosomes named either with or without chr prefix (e.g. 1, chr1)

output will:

- have chromosomes named with chr prefix (e.g. chr1)

Based on process described in https://www.protocols.io/view/genotype-imputation-workflow-v3-0-xbgfijw 

Usage
-----

Turn your microarray data into a block gzipped VCF file with chr prefixed chromosomes and name it `in/input.vcf.gz`

```
docker run \
 --volume `pwd`/in/:/beagle/in/ \
 --volume `pwd`/out/:/beagle/out/ \
 sano-impute:b38-1000genomes-latest \
 /beagle/in/input.vcf.gz \
 /beagle/out/output.vcf.gz
```

Output VCF at `out/output.vcf.gz` and takes about 45 minutes to run on a single sample using a AWS t3a.2xlarge instance (8 vCPU & 32GB).


TODO
----

- X chromosome
- b37 
- b36
- population specifc reference panels
- more reference panel sources

Scripts
-------

`buildall.sh` builds all docker images and tags them as `b38-1000genomes-latest` etc

Docker recap
------------

`docker system df -v` to list current images and sizes

`docker system prune` to cleanup danglers

`docker system prune -a` cleanup everything!