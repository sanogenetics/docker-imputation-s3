Dockerfile for imputation of VCF with beagle

takes an input VCF
assumes that that input VCF has chromosomes named 1,2,3 etc (no "chr" prefix)

bash ref_1_all.sh ref && bash data_1_all.sh in/input.vcf.gz wrk ref  out.vcf.gz