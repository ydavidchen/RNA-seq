#!/bin/bash
# Author: T. Wang, 2017-03
# Maintainer: D. Chen, 2018-04
# This script is to be executed as its entirety (as a `main`) rather than interactively

manual="
Synopsis:
    Composite pipeline (fastqc, STAR/RSEM) for analyzing all RNA-seq FASTQ files in a user-specified directory
Usage:
    executeMe.sh -g <hg19> <-p> -i <path_of_inputs> -o <path_of_outputs> -t <threads>
Options:
    -g    set hg19 as reference genome, default is hg19
    -p    set this for paired-end data
    -i    directory of input fastq or fastq.gz files
    -o    directory for output files, please use full pathway
    -t    average number of threads for each sample, must be integer, default is 1
"

echo "========================== Process begins ========================== "
if [[ $# -le 0 ]]; then
    echo "${manual}" >&2
    exit 1
fi

while getopts :g:pi:o:t: ARGS
do
case $ARGS in
    g)
        genome=$OPTARG
        ;;
    p)
        datatype="PE"
        ;;
    i)
        pathin=$OPTARG
        ;;
    o)
        pathout=$OPTARG
        ;;
    t)
        threads=$OPTARG
        ;;
    :)
        echo "no value for option: $OPTARG"
        echo "${manual}" >&2
        exit 1
        ;;
    *)
        echo "unknow option: $OPTARG"
        echo "${manual}" >&2
        exit 1
        ;;
esac
done

echo "=== Current working directory and settings ==="
cd .
echo $PWD
echo ${genome:="hg19"} >/dev/null
echo ${datatype:="PE"} >/dev/null
echo ${threads:=1} >/dev/null
echo ${pathin:="nopathin"} >/dev/null
echo ${pathout:="nopathout"} >/dev/null
if [[ ${pathin} == "nopathin" ]]; then
    echo "=== ERROR! Please set directory of inputs === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathin} ]]; then
    echo "=== ERROR! Input directory does not exist! === ${manual}" >&2
    exit 1
fi

if [[ ${pathout} == "nopathout" ]]; then
    echo "=== ERROR! Please set directory of outputs === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathout} ]]; then
    echo "=== ERROR! Output directory does not exist! === ${manual}" >&2
    exit 1
fi

(( ${threads} )) 2>/dev/null
if [[ $? != 0 || ! ${threads} -ge 1 ]]; then
    echo "=== ERROR! Thread number must be a positive integer! === ${manual}" >&2
    exit 1
fi
echo ""
echo "=== Set ${threads} threads for each sample in parallel analysis ==="

#------------------------------(Important) Specify genomic annotation paths------------------------------
if [[ ${genome} == "hg19" ]]; then
    echo "=== Set reference genome to: hg19 ==="
    genome_RSEM_STAR_index="/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/RSEM_STAR_Index"
    genome_RSEM_STAR_ref="/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/RSEM_STAR_Index/human_ref"
    genome_fa="/iGenomes/Homo_sapiens/UCSC/hg19/WholeGenomeFasta/genome.fa"
    genome_gtf="/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
    genome_refSeq_Bed="/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_RefSeq.bed"
else
    echo "=== Your reference genome is unsupported! Currently, supported genomes are: hg19 === ${manual}" >&2
    exit 1
fi

#------------------------------Input files check/loading------------------------------
fastq_files=(`find ${pathin} -maxdepth 1 -name "*.fastq*" -type f | sort`)
if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "=== ERROR! No fastq or fastq.gz file found in input directory! === ${manual}" >&2
    exit 1
else
    echo "=== There are ${#fastq_files[@]} raw fastq files in input directory ==="
fi

#------------------------------Key Step: fastQC, STAR Alignment & RSEM Quantification------------------------------
if [[ ${datatype} == "PE" ]]; then
    echo "=== Performing QC with pair-end fastqc... ==="
    mkdir ${pathout}/fastqc
    for file in ${fastq_files[@]}; do fastqc $file -t ${threads} -o ${pathout}/fastqc/ & done
    wait
    fastqc_files=(`ls ${pathout}/fastqc/*fastqc.zip`)
    if [[ ${#fastq_files[@]} -ne ${#fastqc_files[@]} ]]; then
        echo "=== WARNING! Some paired-end samples failed fastqc! ===" >&2
        exit 1
    fi
    multiqc -o ${pathout}/qc/ ${pathout}/fastqc/
    mv ${pathout}/fastqc ${pathout}/qc/
    echo "=== Loading genome into shared memory for parallel alignment... ==="
    mkdir ${pathout}/tmp
    STAR --genomeDir ${genome_RSEM_STAR_index} --genomeLoad LoadAndExit --outFileNamePrefix ${pathout}/tmp/
    rm -rdf ${pathout}/tmp
    echo "=== Parallel analyzing paired-end RNA-seq data...==="
    echo "=== See log.all.txt in each subfolder ==="
    mkdir ${pathout}/process
    for ((i=0; i<${#fastq_files[@]}; i+=2)); do
        subfolder=$(echo ${fastq_files[$i]} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)
        mkdir ${pathout}/process/${subfolder}
        echo "=== Now analyzing sample ${subfolder} ==="
        {
            echo "=== Trimming data and QC for sample ${subfolder} ==="
            trim_galore ${fastq_files[$i]} ${fastq_files[$i+1]} --length 25 --suppress_warn --paired --gzip --o ${pathout}/process/${subfolder}/
            echo "=== Screening contamination for sample ${subfolder} ==="
            fastq_screen ${pathout}/process/${subfolder}/*.fq.gz --subset 1000000 --threads ${threads} --quiet
            echo "=== Mapping to reference genome with STAR for sample ${subfolder} ==="
            STAR --genomeDir ${genome_RSEM_STAR_index} --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --runThreadN ${threads}  --genomeLoad LoadAndKeep  --outSAMtype BAM Unsorted  --quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted  --outFileNamePrefix ${pathout}/process/${subfolder}/  --readFilesCommand zcat  --readFilesIn ${pathout}/process/${subfolder}/*.fq.gz
            echo "=== Assembling transcripts with RSEM based on STAR alignments for sample ${subfolder} ==="
            rsem-calculate-expression --strandedness reverse -p ${threads} -q --time --paired-end --alignments ${pathout}/process/${subfolder}/Aligned.toTranscriptome.out.bam ${genome_RSEM_STAR_ref} ${pathout}/process/${subfolder}/rsem
            echo "=== Sorting genome alignments by coordinate for sample ${subfolder} ==="
            samtools sort -@ ${threads} -m 1G -o ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam ${pathout}/process/${subfolder}/Aligned.out.bam
            samtools index ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam
            echo "=== Annotating bam for sample ${subfolder} ==="
            read_distribution.py -i ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam -r ${genome_refSeq_Bed} > ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam.distribution
            mkdir ${pathout}/process/${subfolder}/tmp_java
            echo "=== Marking duplicates for sample ${subfolder} ==="
            java -Xmx10g -XX:ParallelGCThreads=${threads} -jar $PICARD MarkDuplicates I=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.bam O=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam M=${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=${pathout}/process/${subfolder}/tmp_java/
            echo "=== Getting flag states for sample ${subfolder} ==="
            samtools flagstat ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.bam > ${pathout}/process/${subfolder}/Aligned.sortedByCoord.out.dedupped.flagstat
            rm -rdf ${pathout}/process/${subfolder}/tmp_java
            rm ${pathout}/process/${subfolder}/Aligned.out.bam
        } > ${pathout}/process/${subfolder}/log.all.txt 2>&1 &
    done
wait
else
    echo "=== ERROR! Wrong input of data type (PE only) === ${manual}" >&2
    exit 1
fi

## Implement STAR framework:
echo "=== Analysis complete! Releasing genome from shared memory... ==="
mkdir ${pathout}/tmp
STAR --genomeDir ${genome_RSEM_STAR_index} --genomeLoad Remove --outFileNamePrefix ${pathout}/tmp/
rm -rdf ${pathout}/tmp

#------------------------------Generating & exporting matrices of gene expression results------------------------------
## Generate observed featureCounts:
echo "=== Generating matrix for all samples with featureCounts... ==="
if [[ ${datatype} == "PE" ]]; then
    featureCounts -p -C -s 2 -T 10 -Q 10 -a ${genome_gtf} -o ${pathout}/process/counts.txt ${pathout}/process/*/Aligned.sortedByCoord.out.bam
else
    featureCounts -s 2 -T 10 -Q 10 -a ${genome_gtf} -o ${pathout}/process/counts.txt ${pathout}/process/*/Aligned.sortedByCoord.out.bam
fi
sed -i '1d' ${pathout}/process/counts.txt
cut -f 1,7- ${pathout}/process/counts.txt > ${pathout}/process/counts2.txt
mv ${pathout}/process/counts2.txt ${pathout}/process/counts.txt

## Generate matrix of expected counts:
echo "=== Generating featureCount matrix... ==="
rsem-generate-data-matrix ${pathout}/process/*/rsem.genes.results > ${pathout}/process/rsem_expected_count_genes.txt
rsem-generate-data-matrix ${pathout}/process/*/rsem.isoforms.results > ${pathout}/process/rsem_expected_count_isoforms.txt

## Generate FPKM matrix:
echo "=== Generating FPKM matrix for all samples... ==="
folders=(`ls -d ${pathout}/process/*/ | rev | cut -d '/' -f2 | rev`)
mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,7 ${pathout}/process/${folders[$i]}/rsem.genes.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/fpkm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/fpkm_all.txt > ${pathout}/process/tmp/fpkm_table.txt
headarr=('gene_id' 'transcript_id(s)' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/fpkm_table.txt > ${pathout}/process/rsem_fpkm_genes.txt
rm -rdf ${pathout}/process/tmp

mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,7 ${pathout}/process/${folders[$i]}/rsem.isoforms.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/fpkm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/fpkm_all.txt > ${pathout}/process/tmp/fpkm_table.txt
headarr=('transcript_id' 'gene_id' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/fpkm_table.txt > ${pathout}/process/rsem_fpkm_isoforms.txt
rm -rdf ${pathout}/process/tmp

## Generate TPM matrix: 
echo "=== Generating TPM matrix for all samples... ==="
folders=(`ls -d ${pathout}/process/*/ | rev | cut -d '/' -f2 | rev`)
mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,6 ${pathout}/process/${folders[$i]}/rsem.genes.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/tpm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/tpm_all.txt > ${pathout}/process/tmp/tpm_table.txt
headarr=('gene_id' 'transcript_id(s)' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/tpm_table.txt > ${pathout}/process/rsem_tpm_genes.txt
rm -rdf ${pathout}/process/tmp

mkdir ${pathout}/process/tmp
for ((i=0; i<${#folders[@]}; i++)); do cut -f 1,2,6 ${pathout}/process/${folders[$i]}/rsem.isoforms.results > ${pathout}/process/tmp/${folders[$i]}; done
for file in ${pathout}/process/tmp/*; do sed -i '1d' $file; done
for file in ${pathout}/process/tmp/*; do sort -k1,1 -k2,2 $file > $file.sort; done
paste ${pathout}/process/tmp/*.sort > ${pathout}/process/tmp/tpm_all.txt
column=$(seq 3 3 $((3*${#folders[@]})) | paste -sd ',')
cut -f 1-${column} ${pathout}/process/tmp/tpm_all.txt > ${pathout}/process/tmp/tpm_table.txt
headarr=('transcript_id' 'gene_id' ${folders[*]})
printf $'%s\t' ${headarr[@]} | cut -f 1-${#headarr[@]} > ${pathout}/process/tmp/header.txt
cat ${pathout}/process/tmp/header.txt ${pathout}/process/tmp/tpm_table.txt > ${pathout}/process/rsem_tpm_isoforms.txt
rm -rdf ${pathout}/process/tmp

## Generate summary table:
echo "=== Generating summary table... ==="
summarize_RNAseq_STAR.RSEM.pl ${pathout}/process/ ${pathout}/process/summary.txt # summarize_RNAseq_RSEM.pl is sudo-cp'ed to /usr/local/bin/
mkdir ${pathout}/tables
mv ${pathout}/process/*txt* ${pathout}/tables/

echo "========================== Process finished ========================== "
