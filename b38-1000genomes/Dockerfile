# openJDK on Debial
FROM openjdk:19-bullseye

RUN mkdir -p /beagle/ref && mkdir -p /beagle/bin
WORKDIR /beagle

COPY ref_10_deps.sh /beagle/
RUN bash /beagle/ref_10_deps.sh

COPY ref_15_jars.sh /beagle/
RUN bash /beagle/ref_15_jars.sh

COPY ref_20_map.sh /beagle/
RUN bash /beagle/ref_20_map.sh

# this will take a while and use a lot of space and bandwidth
COPY ref_30_seqs.sh /beagle/
COPY ref_31_seq.sh /beagle/
RUN bash /beagle/ref_30_seqs.sh

# this will take a while and use a lot of space and bandwidth
COPY ref_40_chroms.sh /beagle/
COPY ref_41_chrom.sh /beagle/
RUN bash /beagle/ref_40_chroms.sh

# add the scripts that will be used to run the imputation
COPY data_1_all.sh /beagle/
COPY data_2_chr.sh /beagle/

ENTRYPOINT ["/beagle/data_1_all.sh"]
