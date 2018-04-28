######################################################################################################
# Generate genomeParameter.txt file
# Script author: David Chen (@ydavidchen)
# Date: April 28, 2018
# Notes:
# 1. This script is to be run interactively in MacOS Terminal. 
######################################################################################################

## Check if files exist:
if [ -f /Users/DavidKevinChen/repos/RNA-seq/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_RefSeq.bed ] #absolute path required
# /Users/DavidKevinChen/repos/RNA-seq/iGenomes/Homo_sapiens/UCSC/hg19/WholeGenomeFasta/genome.fa #also test
# /Users/DavidKevinChen/repos/RNA-seq/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf #also test
then
    echo "The file exists. OK to proceed."
else
    echo "The file does not exist or has wrong path!"
fi

## Assume you're in the RNA-seq repo
cd iGenomes/Homo_sapiens/UCSC/hg19/

mkdir GENOME_data #grants permission also
mkdir GENOME_data/star #grants permission also

## Prepare genomeParameter.txt etc:
STAR --runThreadN 40 \
 --runMode genomeGenerate \
 --genomeDir GENOME_data/star \
 --genomeFastaFiles WholeGenomeFasta/genome.fa \
 --sjdbGTFfile Annotation/Genes/genes.gtf \
 --limitGenomeGenerateRAM=85543555797 #can try 124544990592

 ## Copy the files to required directory:
cp GENOME_data/star/genomeParameters.txt Sequence/RSEM_STAR_Index/

