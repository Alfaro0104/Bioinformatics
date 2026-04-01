#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o assembly.test.log
#SBATCH --account=dtalfaro5077
#SBATCH --partition=silver

#load tools
module load biological/samtools_1.23
module load biological/java
module load biological/perl_5.40


#sets working directory and sample name
export PROJ_DIR=/export/home/bio_class/dtalfaro5077/Bioinformatics_Lab_9
cd $PROJ_DIR
export SRR=SRR5324768

#makes needed folders
mkdir -p alignment
mkdir -p genome
mkdir -p variants

#makes sequence dictionary if it does not already exist
if [ ! -f genome/Thermus_thermophilus_TTHNAR1.dict ]; then
    java -jar /export/share/software/biological/picard/picard.jar \
        CreateSequenceDictionary \
        REFERENCE=genome/Thermus_thermophilus_TTHNAR1.fa \
        OUTPUT=genome/Thermus_thermophilus_TTHNAR1.dict
fi

#builds the bowtie2 index
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2-build \
    genome/Thermus_thermophilus_TTHNAR1.fa \
    genome/Thermus_thermophilus_TTHNAR1

#makes the SAM file
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 \
    -x genome/Thermus_thermophilus_TTHNAR1 \
    -1 fastq/${SRR}_pass_1.fastq.gz \
    -2 fastq/${SRR}_pass_2.fastq.gz \
    --sensitive-local \
    --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
    > alignment/${SRR}.sam

#converts the SAM to BAM
samtools view -hb alignment/${SRR}.sam \
    | samtools sort -l 5 -o alignment/${SRR}.bam

#makes pileup from BAM
samtools mpileup -f ncbi_dataset/data/GCA_900604845.1/GCA_900604845.1_TTHNAR1_genomic.fna \
    alignment/${SRR}.bam > variants/${SRR}.pileup

#makes the VCF from pileup
echo "##fileformat=VCFv4.2" > variants/${SRR}.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> variants/${SRR}.vcf
head -n 20 variants/${SRR}.pileup | awk '{print $1"\t"$2"\t.\t"$3"\tN\t.\t.\t."}' >> variants/${SRR}.vcf

#makes a consensus file instead of a vcf
samtools consensus -f fasta -o ${SRR}_consensus.fasta alignment/${SRR}.bam
