#!/bin/bash
# This script uses the SRA toolkit (NCBI, SRA Toolkit Development Team | https://github.com/ncbi/sra-tools) to download the raw fastq data from the studies
####################################################################

# command line input of accession list and output path:
AccList=$1
out=$2

# prefetching:
for SRR in $(cat ${AccList}); do
        mkdir -p ${out}/${SRR}
        prefetch ${SRR} -O ${out}
done
echo "prefetching done!"

# copying accession list to SRR directory; can't do before bc prefetch needs to be pointed to an empty folder.
cp ${AccList} ${out}

# getting the fastq files:
cd ${out}
for SRR in $(cat ${AccList}); do
        fasterq-dump ${SRR}

done
echo "loading done!"