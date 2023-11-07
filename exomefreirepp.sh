#instalar wget
brew instal wget

#REFERENCE GENOME (download, index, create dictionary)

#download the reference genome
wget [URL_of_the_reference_genome]

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#to index the reference genome 

bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#create a dictionary of reference genome
java -jar /home/freirepp/software/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.primary_assembly.fa O=Homo_sapiens.GRCh38.dna.primary_assembly.dict

#after indexing and creating the dictionary, the following files should appear in your directory:

-rw-rw-r-- 1 freirepp freirepp  31K set 29 16:03 Homo_sapiens.GRCh38.dna.primary_assembly.dict
-rw-rw-r-- 1 freirepp freirepp 3,0G abr 21 13:28 Homo_sapiens.GRCh38.dna.primary_assembly.fa
-rw-rw-r-- 1 freirepp freirepp 6,3K set 29 15:22 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
-rw-rw-r-- 1 freirepp freirepp  18K ago 18 15:26 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.amb
-rw-rw-r-- 1 freirepp freirepp  17K ago 18 15:26 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.ann
-rw-rw-r-- 1 freirepp freirepp 2,9G ago 18 15:25 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.bwt
-rw-rw-r-- 1 freirepp freirepp 6,3K set 29 13:06 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai
-rw-rw-r-- 1 freirepp freirepp 740M ago 18 15:26 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.pac
-rw-rw-r-- 1 freirepp freirepp 1,5G ago 18 15:45 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.sa

#to install FASTQC
# nohup means ...(para rodar no servidor sem vc está na sua conta)
# & means ... (& no final do nohup, para o nohup funcionar)
nohup conda install -c bioconda fastqc &

#to remove the low quality sequences (I used cutadapt, but we have the option of trimmomatic)
#to install trimmomatic 
conda install -c bioconda trimmomatic

#to analyze the sequence quality using FASTQC and create a html file to analyze the quality of your sequence
#here, -o will create a directory containing the results;
# -t 3--nogroup ???; 
# *.fastq.gz: all the files that finalize with fastq.gz
fastqc -o ./fastqc_before -t 3 --nogroup *.fastq.gz

#to install PICARD 
conda install -c bioconda picard
#tive problemas em usar o picard via anaconda, então baixei direto do link do site do picard (https://github.com/broadinstitute/picard/releases/tag/3.1.0)
wget https://github.com/broadinstitute/picard/releases/tag/3.1.0#:~:text=3-,picard.jar,-59.6%20MB

chmod 777 picard-1.119.jar

#Problem: JAVA 1.8 é requerido  
#atualizar java 
conda install openjdk=8

#to align your sequences with the indexing genome [BWA MEM]
bwa mem /home/freirepp/genome_reference/genome_fa/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz *R1* *R2* > BWA_P314_mem.sam

#To check if your SAM file is correct [PICARD]:
java -jar /home/freirepp/software/picard.jar ValidateSamFile I=BWA_P314_mem.sam MODE=SUMMARY
#ERROR:MISSING_READ_GROUP

#Add READER GROUP in SAM file (output do BWA) [PICARD]
java -jar ~/software/picard.jar AddOrReplaceReadGroups  I=BWA_P314_mem.sam O=P314_out_rg.sam RGID=group1 RGLB=library1 RGPL=illumina RGPU=unit1 RGSM=sample1  CREATE_INDEX=true

#To check if your SAM file is correct [PICARD]
java -jar /home/freirepp/software/picard.jar ValidateSamFile I=P314_out_rg.sam MODE=SUMMARY

#to transform SAM file in BAM file [SAMTOOLS]
#you need to be inside the directory that you put your SAM file (ex. Samples_HF)
/home/freirepp/software/samtools-1.18/samtools view -H BWA_P314_mem.sam 
#problem: [main_samview] fail to read the header from "BWA_P314_mem.sam".
#we tried to align again
bwa mem /home/freirepp/genome_reference/genome_fa/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz *R1* *R2* > BWA_P314_mem.sam

#than: 
/home/freirepp/software/samtools-1.18/samtools view -h P314_out_rg.sam -o P314_out_rg.bam 

#to transform .BAM file in .BAI file (it means like to index or mapping the BAM file)[SAMTOOLS]
/home/freirepp/software/samtools-1.18/samtools index BWA_P314_mem.bam

#in your directory you will find the files: 
-rw-rw-r-- 1 freirepp freirepp       2164 set 15 10:40 build.sh
-rw-rw-r-- 1 freirepp freirepp          0 set 15 13:21 BWA_P314.bam
-rw-rw-r-- 1 freirepp freirepp       2841 set 21 17:11 BWA_P314_mem.bam
-rw-rw-r-- 1 freirepp freirepp       1568 set 22 13:43 BWA_P314_mem.bam.bai
-rw-rw-r-- 1 freirepp freirepp 1309745519 set 21 16:04 BWA_P314_mem.sam
drwxrwxr-x 2 freirepp freirepp        209 set 14 14:50 fastqc_before
-rw-rw-r-- 1 freirepp freirepp  161321728 set 13 16:54 P314_HFBP031_S3_L001_R1_001.fastq.gz
-rw-rw-r-- 1 freirepp freirepp  167961960 set 13 16:54 P314_HFBP031_S3_L001_R2_001.fastq.gz

# Create a Sorted Header [PICARD]
# Exception in thread "main" picard.PicardException: This program requires input that are either coordinate or query sorted (according to the header, or at least ASSUME_SORT_ORDER and the content.) 
#Found ASSUME_SORT_ORDER=null and header sortorder=unsorted

/home/freirepp/software/samtools-1.18/samtools sort -o P314_out_rg_sorted.bam P314_out_rg.bam

# Mark duplicates with PICARD
java -jar /home/freirepp/software/picard.jar MarkDuplicates -I P314_out_rg_sorted.bam -O P314_out_dup_rg.bam -M marked_dup_metrics.txt

#indexar o BAM [SAMTOOLS]
/home/freirepp/software/samtools-1.18/samtools index P314_out_dup_rg.bam

# Remove duplicates with PICARD
java -jar /home/freirepp/software/picard.jar MarkDuplicates -I P314_out_rg_sorted.bam -O P314_out_dup_rg.bam -M marked_dup_metrics.txt --REMOVE_DUPLICATES true

# Call variants with GATK - To generate de VCF file with GATK [GATK]

~/software/gatk-4.4.0.0/gatk HaplotypeCaller -R /home/freirepp/genome_reference/genome_fa/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I P314_output.no.duplicates.bam -O P314_variants.vcf

#to transform em VCF [PICARD]
java -jar ~/software/picard.jar AddOrReplaceReadGroups  I=P314_output_no_duplicates.bam O=P314_output_no_duplicates.bam_rg.bam RGID=group1 RGLB=library1 RGPL=illumina RGPU=unit1 RGSM=sample1  CREATE_INDEX=true

#to annotate the variants [ANNOVAR]
upload the VCF file in the website: https://wannovar.wglab.org/

#to convert BAM file to BED file [BEDTOOLS]
bedtools bamtobed -i P314_out_dup_rg.bam > P314.bed
bedtools bamtobed -i P314_out_dup_rg.bam -ed > P314_score.bed 


#################################################

# OPTION: Mark duplicates with Samtools [SAMTOOLS]
/home/freirepp/software/samtools-1.18/samtools markdup -r BWA_P314_mem.bam marked_output.bam

#create a new envviroment (envs)
conda activate java17

#para entrar no environment java17
conda activate java17

#install new version of JAVA:
conda install -c conda-forge openjdk

* (tudo) ex.: *fastq.gz (Tudo que terminar com fastq.gz)






