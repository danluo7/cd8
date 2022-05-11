# SK2 knockout

sk2 knockout vs. wildtype


# raw data processing


## check quality of reads

    fastqc *.fastq.gz 

    python3 -m multiqc . 

    mkdir fastqc
    mv *fastqc* fastqc


For mouse mm10 genome hisat2 pre-built index:

    mkdir -p cd8/RNA_REF_FA

set environmental variables by nano .bashrc at the home directory, then setting cd8=/home/floyd/ubuntu/workspace/cd8, and setting cd8_data=/media/floyd/OneTouch/SK2RNAseq


    cd $cd8/RNA_REF_FA
    wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
    tar -xzvf mm10_genome.tar.gz


download reference genome

    mkdir -p $cd8/RNA_REF_GTF
    cd $cd8/RNA_REF_GTF
    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
    gzip -d mm10.ncbiRefSeq.gtf.gz
    
    
## alignment using hisat2

take a look at the files to see the platform info to find the barcode index:

        gzip -cd Wild-CDT1_S1_L001_R1_001.fastq.gz | head
        gzip -cd Wild-CDT2_S2_L001_R1_001.fastq.gz | head
        gzip -cd Wild-CDT3_S3_L001_R1_001.fastq.gz | head
       
        gzip -cd Sk2KOCD1_S4_L001_R1_001.fastq.gz | head
        gzip -cd Sk2KOCD2_S5_L001_R1_001.fastq.gz | head
        gzip -cd Sk2KOCD3_S6_L001_R1_001.fastq.gz | head
        
        gzip -cd MDSCsWild1_S7_L001_R1_001.fastq.gz | head
        gzip -cd MDSCsWild2_S8_L001_R1_001.fastq.gz | head
        gzip -cd MDSCsWild3_S9_L001_R1_001.fastq.gz | head
       
        gzip -cd MDSCsSK2KO1_S10_L001_R1_001.fastq.gz | head
        gzip -cd MDSCsSK2KO2_S11_L001_R1_001.fastq.gz | head
        gzip -cd MDSCsSK2KO3_S12_L001_R1_001.fastq.gz | head
        

output: 
@A00257:731:HV577DRXY:1:2101:1253:1000 1:N:0:TGTCGCTGGT+AGTCAGACGT
NCACTCTTTGGTTCCAGGAAACCCCGGCTCCCAATCAGCCCCGTGTGCTTC

    mkdir -p $cd8_data/alignments

### At this point decide on naming conventions. Order the samples logically

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_TATGCCTTAC+TAATGTGTCT --rg PL:ILLUMINA --rg PU:HV577DRXY.1.TATGCCTTAC+TAATGTGTCT -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/Wild-CDT1_S1_L001_R1_001.fastq.gz -2 $cd8_data/Wild-CDT1_S1_L001_R2_001.fastq.gz -S $cd8_data/alignments/1_WTCD8.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_ACCGTTACAA+TCGCTATGAG --rg PL:ILLUMINA --rg PU:HV577DRXY.1.ACCGTTACAA+TCGCTATGAG -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/Wild-CDT2_S2_L001_R1_001.fastq.gz -2 $cd8_data/Wild-CDT2_S2_L001_R2_001.fastq.gz -S $cd8_data/alignments/2_WTCD8.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_TATGCCTTAC+TAATGTGTCT --rg PL:ILLUMINA --rg PU:HV577DRXY.1.TATGCCTTAC+TAATGTGTCT -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/Wild-CDT3_S3_L001_R1_001.fastq.gz -2 $cd8_data/Wild-CDT3_S3_L001_R2_001.fastq.gz -S $cd8_data/alignments/3_WTCD8.sam


    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_ACAAGTGGAC+AACATCGCGC --rg PL:ILLUMINA --rg PU:HV577DRXY.1.ACAAGTGGAC+AACATCGCGC -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/Sk2KOCD1_S4_L001_R1_001.fastq.gz -2 $cd8_data/Sk2KOCD1_S4_L001_R2_001.fastq.gz -S $cd8_data/alignments/4_SK2KOCD8.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_TGGTACCTAA+AGTACTCATG --rg PL:ILLUMINA --rg PU:HV577DRXY.1.TGGTACCTAA+AGTACTCATG -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/Sk2KOCD2_S5_L001_R1_001.fastq.gz -2 $cd8_data/Sk2KOCD2_S5_L001_R2_001.fastq.gz -S $cd8_data/alignments/5_SK2KOCD8.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_TTGGAATTCC+GTATTGACGT --rg PL:ILLUMINA --rg PU:HV577DRXY.1.TTGGAATTCC+GTATTGACGT -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/Sk2KOCD3_S6_L001_R1_001.fastq.gz -2 $cd8_data/Sk2KOCD3_S6_L001_R2_001.fastq.gz -S $cd8_data/alignments/6_SK2KOCD8.sam


    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_CCTCTACATG+AGGAGGTATC --rg PL:ILLUMINA --rg PU:HV577DRXY.1.CCTCTACATG+AGGAGGTATC -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/MDSCsWild1_S7_L001_R1_001.fastq.gz -2 $cd8_data/MDSCsWild1_S7_L001_R2_001.fastq.gz -S $cd8_data/alignments/7_MDSCWT.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_GGAGCGTGTA+ACTTACGGAT --rg PL:ILLUMINA --rg PU:HV577DRXY.1.GGAGCGTGTA+ACTTACGGAT -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/MDSCsWild2_S8_L001_R1_001.fastq.gz -2 $cd8_data/MDSCsWild2_S8_L001_R2_001.fastq.gz -S $cd8_data/alignments/8_MDSCWT.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_GTCCGTAAGC+AAGATACACG --rg PL:ILLUMINA --rg PU:HV577DRXY.1.GTCCGTAAGC+AAGATACACG -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/MDSCsWild3_S9_L001_R1_001.fastq.gz -2 $cd8_data/MDSCsWild3_S9_L001_R2_001.fastq.gz -S $cd8_data/alignments/9_MDSCWT.sam


    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_ACTTCAAGCG+TTCATGGTTC --rg PL:ILLUMINA --rg PU:HV577DRXY.1.ACTTCAAGCG+TTCATGGTTC -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/MDSCsSK2KO1_S10_L001_R1_001.fastq.gz -2 $cd8_data/MDSCsSK2KO1_S10_L001_R2_001.fastq.gz -S $cd8_data/alignments/10_MDSCSK2KO.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_TCAGAAGGCG+TATGATGGCC --rg PL:ILLUMINA --rg PU:HV577DRXY.1.TCAGAAGGCG+TATGATGGCC -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/MDSCsSK2KO2_S11_L001_R1_001.fastq.gz -2 $cd8_data/MDSCsSK2KO2_S11_L001_R2_001.fastq.gz -S $cd8_data/alignments/11_MDSCSK2KO.sam

    hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_GCGTTGGTAT+GGAAGTATGT --rg PL:ILLUMINA --rg PU:HV577DRXY.1.GCGTTGGTAT+GGAAGTATGT -x $cd8/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $cd8_data/MDSCsSK2KO3_S12_L001_R1_001.fastq.gz -2 $cd8_data/MDSCsSK2KO3_S12_L001_R2_001.fastq.gz -S $cd8_data/alignments/12_MDSCSK2KO.sam


## Convert SAM files into BAM files

    cd $cd8_data/alignments
    
    samtools sort -@ 10 -o 1_WTCD8.bam 1_WTCD8.sam
    samtools sort -@ 10 -o 2_WTCD8.bam 2_WTCD8.sam
    samtools sort -@ 10 -o 3_WTCD8.bam 3_WTCD8.sam
   
    samtools sort -@ 10 -o 4_SK2KOCD8.bam 4_SK2KOCD8.sam 
    samtools sort -@ 10 -o 5_SK2KOCD8.bam 5_SK2KOCD8.sam
    samtools sort -@ 10 -o 6_SK2KOCD8.bam 6_SK2KOCD8.sam
    
    samtools sort -@ 10 -o 7_MDSCWT.bam 7_MDSCWT.sam
    samtools sort -@ 10 -o 8_MDSCWT.bam 8_MDSCWT.sam
    samtools sort -@ 10 -o 9_MDSCWT.bam 9_MDSCWT.sam
    
    samtools sort -@ 10 -o 10_MDSCSK2KO.bam 10_MDSCSK2KO.sam
    samtools sort -@ 10 -o 11_MDSCSK2KO.bam 11_MDSCSK2KO.sam
    samtools sort -@ 10 -o 12_MDSCSK2KO.bam 12_MDSCSK2KO.sam
    
    
## Index all bam files for IGV visualization

Make sure that the only files in the directory are the sam and bam/bai files, then:

can do: samtools index 1.bam

or:

    find *.bam -exec echo samtools index {} \; | sh

## use samtoools flagstat to geta basic sumary of an alignment

    mkdir flagstat
    samtools flagstat 1_WTCD8.bam > flagstat/1.bam.flagstat
    samtools flagstat 2_WTCD8.bam > flagstat/2.bam.flagstat

viewing them:

    cat flagstat/1.bam.flagstat


## Use Stringtie to generate expression estimates (FPKM and TPM) from the SAM/BAM files generated by HISAT2 and samtools


    cd /home/floyd/ubuntu/bin
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz
    tar -xzvf stringtie-2.2.1.Linux_x86_64.tar.gz

    ln -s /home/floyd/ubuntu/bin/stringtie-2.1.4.Linux_x86_64/stringtie
    
    mkdir -p $cd8/expression/stringtie/ref_only/
    cd $cd8/expression/stringtie/ref_only/
    
    
stringtie did not install on this installation of ubuntu (had issues with linking/putting stringtie in the gernal bin folder. Chose to run stringtie in the local folder, then just copy the resulting counts files over to the right directory later. 

    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 1_WTCD8/transcripts.gtf -A 1_WTCD8/gene_abundances.tsv $cd8_data/alignments/1_WTCD8.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 2_WTCD8/transcripts.gtf -A 2_WTCD8/gene_abundances.tsv $cd8_data/alignments/2_WTCD8.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 3_WTCD8/transcripts.gtf -A 3_WTCD8/gene_abundances.tsv $cd8_data/alignments/3_WTCD8.bam

    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 4_SK2KOCD8/transcripts.gtf -A 4_SK2KOCD8/gene_abundances.tsv $cd8_data/alignments/4_SK2KOCD8.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 5_SK2KOCD8/transcripts.gtf -A 5_SK2KOCD8/gene_abundances.tsv $cd8_data/alignments/5_SK2KOCD8.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 6_SK2KOCD8/transcripts.gtf -A 6_SK2KOCD8/gene_abundances.tsv $cd8_data/alignments/6_SK2KOCD8.bam


    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 7_MDSCWT/transcripts.gtf -A 7_MDSCWT/gene_abundances.tsv $cd8_data/alignments/7_MDSCWT.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 8_MDSCWT/transcripts.gtf -A 8_MDSCWT/gene_abundances.tsv $cd8_data/alignments/8_MDSCWT.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 9_MDSCWT/transcripts.gtf -A 9_MDSCWT/gene_abundances.tsv $cd8_data/alignments/9_MDSCWT.bam

    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 10_MDSCSK2KO/transcripts.gtf -A 10_MDSCSK2KO/gene_abundances.tsv $cd8_data/alignments/10_MDSCSK2KO.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 11_MDSCSK2KO/transcripts.gtf -A 11_MDSCSK2KO/gene_abundances.tsv $cd8_data/alignments/11_MDSCSK2KO.bam
    ./stringtie --fr -p 10 -G $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 12_MDSCSK2KO/transcripts.gtf -A 12_MDSCSK2KO/gene_abundances.tsv $cd8_data/alignments/12_MDSCSK2KO.bam


## htseq for raw counts DE

    mkdir -p $cd8/expression/htseq_counts
    cd $cd8/expression/htseq_counts

    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/1_WTCD8.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 1_WTCD8.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/2_WTCD8.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 2_WTCD8.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/3_WTCD8.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 3_WTCD8.tsv

    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/4_Sk2KOCD8.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 4_Sk2KOCD8.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/5_Sk2KOCD8.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 5_Sk2KOCD8.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/6_Sk2KOCD8.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 6_Sk2KOCD8.tsv

    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/7_MDSCWT.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 7_MDSCWT.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/8_MDSCWT.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 8_MDSCWT.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/9_MDSCWT.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 9_MDSCWT.tsv

    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/10_MDSCSK2KO.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 10_MDSCSK2KO.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/11_MDSCSK2KO.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 11_MDSCSK2KO.tsv
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $cd8_data/alignments/12_MDSCSK2KO.bam $cd8/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 12_MDSCSK2KO.tsv
    
    
    
## Use Ballgown in R for differential expression (DE) analysis (then PCA) using output from Stringtie

    mkdir -p $cd8/de/ballgown/ref_only
	cd $cd8/de/ballgown/ref_only/


Use printf to create/print a table with ids, type (each type of sample is a type), and path to the file, as the header. Then n returns a new line. 
## (note: ids needs to match the file folder names created by stringtie)

Bascially, need a table that needs to look like this to feed into R:

ids type path-to-file-011_invitro_1 011 $gbm/expression/stringtie/1 011_invitro_2 011 $gbm/expression/stringtie/2 ... ...

goal is to generate a header file to load into R, for ALL samples for principal component analysis (the simplest form of multidimentional scaling), and also a file for pairwise comparisons. since we have a ton of comparisisons, might just not do this for now and only do the PCA. 

file for all 011 samples for PCA: (this is how the script should look like (without the enters inbetween each line):

printf "\"ids\",\"type\",\"path
\"\n\"1_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/1_WTCD8
\"\n\"2_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/2_WTCD8
\"\n\"3_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/3_WTCD8

\"\n\"4_SK2KOCD8\",\"SK2KO_CD8\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8
\"\n\"5_SK2KOCD8\",\"SK2KO_CD8\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8
\"\n\"6_SK2KOCD8\",\"SK2KO_CD8\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8

\"\n\"7_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT
\"\n\"8_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT
\"\n\"9_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT

\"\n\"10_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO
\"\n\"11_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO
\"\n\"12_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO
\"\n" > SK2_all.csv

	printf "\"ids\",\"type\",\"path\"\n\"1_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/1_WTCD8\"\n\"2_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/2_WTCD8\"\n\"3_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/3_WTCD8\"\n\"4_SK2KOCD8\",\"SK2KO_CD8\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8\"\n\"5_SK2KOCD8\",\"SK2KO_CD8\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8\"\n\"6_SK2KOCD8\",\"SK2KO_CD8\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8\"\n\"7_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT\"\n\"8_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT\"\n\"9_MDSCWT\",\"MDSC_WT\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT\"\n\"10_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO\"\n\"11_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO\"\n\"12_MDSCSK2KO\",\"MDSC_SK2KO\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO\"\n" > SK2_all.csv


R script:


	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("SK2_all.csv")  


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg


	bg_table = texpr(bg, 'all')


	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)


	save(bg, file='bg.rda')

	bg
	
	
	
	pdf(file="SK2_R_output.pdf")

	working_dir = "~/workspace/cd8/de/ballgown/ref_only"
	setwd(working_dir)
	dir()


Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
	
	load('bg.rda')


Load gene names for lookup later in the tutorial
	
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])


Pull the gene_expression data frame from the ballgown object

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

View the column names

	colnames(gene_expression)

View the row names
	
	head(row.names(gene_expression))

Determine the dimensions of the dataframe.  'dim()' will return the number of rows and columns

	dim(gene_expression)


Just for fun, check BRD4 expression across all 8 samples:

	i = row.names(gene_expression) == "Brd4"
	gene_expression[i,]



Load the transcript to gene index from the ballgown object. Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes  


	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
	
	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count


> length(row.names(transcript_gene_table)) #Transcript count
[1] 106520
> length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count
[1] 35976
> 


Plot the number of transcripts per gene. Many genes will have only 1 transcript, some genes will have several transcripts. Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol). Then use the 'hist' command to create a histogram of these counts

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)


Plot the distribution of transcript sizes as a histogram. lengths will be those of known transcripts. Good QC step: we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts.

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

View the summary FPKM values (minimum and maximum FPKM values) for any particular library

	min(gene_expression[,"FPKM.1_WTCD8"])
	max(gene_expression[,"FPKM.3_WTCD8"])


Set the minimum non-zero FPKM values by one of two ways:

coverting 0's to NA, and calculating the minimum or all non NA values
two ways: 
zz = fpkm_matrix[,data_columns]
zz[zz==0] = NA
min_nonzero = min(zz, na.rm=TRUE)
min_nonzero


Alternatively just set min value to 1
	
	min_nonzero=1

Set the columns for finding FPKM and create shorter names for figures

	data_columns=c(1:12)
	short_names=c("WT_CD8_1","WT_CD8_2","WT_CD8_3","SK2KO_CD8_1","SK2KO_CD8_2","SK2KO_CD8_3","MDSC_WT_1","MDSC_WT_2","MDSC_WT_3","MDSC_SK2KO_1","MDSC_SK2KO_2","MDSC_SK2KO_3")


Plot range of values and general distribution of FPKM values for all 8 libraries

Create boxplots using different colors by setting storing the colors of the columns in a variable called data_colors. then display on a log2 scale and add the minimum non-zero value to avoid log2(0). Note that the bold horizontal line on each boxplot is the median.


	colors()
	data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3","seagreen1","seagreen2","seagreen3","red1","red2","red3")

	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 12 sample libraries")


## Compare the correlation distance between all replicates

Calculate the FPKM sum for all 12 libraries

	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

Filter out genes with a grand sum FPKM of less than 12

	i = which(gene_expression[,"sum"] > 10)


Calculate the correlation between all pairs of data

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r
	
	
## Plot MDS.
Convert correlation to distance, and use 'multi-dimensional scaling' to plot the relative differences between libraries, by calculating 2-dimensional coordinates to plot points for each library using eigenvectors (eig=TRUE). d, k=2 means 2 dimensions
	
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.004,0.004), ylim=c(-0.0015,0.0015))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))



close out the PDF
dev.off()



## comparing CD8 to non-CD8, using stattest function, and output a heatmap

#Need to make a new header file and recode the "type" header,since this will determine what gets compared in a stattest. the stattest function will use type as a covariate and use fpkm as a meansurement. since this function can't compare multiple things, need to make another file called GBM049_all_stattest.csv and make the "type" tissue vs. non-tissue. Then heatmaps can be done comparing invitro to everything else.


	mkdir -p $cd8/de/ballgown/ref_only/stattest
	cd $cd8/de/ballgown/ref_only/stattest
	
printf "\"ids\",\"type\",\"path
\"\n\"1_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/1_WTCD8
\"\n\"2_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/2_WTCD8
\"\n\"3_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/3_WTCD8

\"\n\"4_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8
\"\n\"5_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8
\"\n\"6_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8

\"\n\"7_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT
\"\n\"8_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT
\"\n\"9_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT

\"\n\"10_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO
\"\n\"11_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO
\"\n\"12_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO


	printf "\"ids\",\"type\",\"path\"\n\"1_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/1_WTCD8\"\n\"2_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/2_WTCD8\"\n\"3_WTCD8\",\"WT_CD8\",\"$cd8/expression/stringtie/ref_only/3_WTCD8\"\n\"4_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/4_SK2KOCD8\"\n\"5_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/5_SK2KOCD8\"\n\"6_SK2KOCD8\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/6_SK2KOCD8\"\n\"7_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/7_MDSCWT\"\n\"8_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/8_MDSCWT\"\n\"9_MDSCWT\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/9_MDSCWT\"\n\"10_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/10_MDSCSK2KO\"\n\"11_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/11_MDSCSK2KO\"\n\"12_MDSCSK2KO\",\"non_WT_CD8\",\"$cd8/expression/stringtie/ref_only/12_MDSCSK2KO\"\n" > SK2_all_stattest.csv
	
	
#now rerun all the R scripts to do stattest and heat map.



	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("SK2_all_stattest.csv")  

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')

	bg

	pdf(file="SK2_R_output_stattest.pdf")


	dir()

	load('bg.rda')

	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

	colnames(gene_expression)
	dim(gene_expression)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]

	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)

	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

	min(gene_expression[,"FPKM.1_WTCD8"])
	max(gene_expression[,"FPKM.2_WTCD8"])

	min_nonzero=1
	
	
#make sure the columns AND short names are correct, this is for the names of columns in the heatmap 


	data_columns=c(1:12) 		
	short_names=c("CD8_WT","CD8_WT","CD8_WT","CD8_SK2KO","CD8_SK2KO","CD8_SK2KO","MDSC_WT","MDSC_WT","MDSC_WT","MDSC_SK2KO","MDSC_SK2KO","MDSC_SK2KO")


#colors()

data_colors=c("black","black","black","green","green","green","red","red","red","blue","blue","blue")


boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 12 samples")





gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 10)

r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r



## Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


## View the distribution of differential expression values as a histogram

#Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) CD8_WT vs non-CD8WT", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)
	
	#Display the grand expression values from WTCD8 vs non-WTCD8 and mark those that are significantly differentially expressed. Make sure all the columns are correctly indicated. 

	gene_expression[,"CD8_WT"]=apply(gene_expression[,c(1:3)], 1, mean)
	gene_expression[,"non-CD8_WT"]=apply(gene_expression[,c(4:12)], 1, mean)

	x=log2(gene_expression[,"CD8_WT"]+min_nonzero)
	y=log2(gene_expression[,"non-CD8_WT"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="FPKM (log2)", ylab="FPKM (log2)", main="tissue vs invitro FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



Get the gene symbols for the top 25 genes (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown_CD8_WT_vs_non-CD8_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output

	output[1:25,c(1,4,5)]
	
	
	
	
## Create a heatmap to vizualize expression differences between the 16 samples
	

#Define custom dist and hclust functions for use with heatmaps

	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="average")}

	main_title="sig DE Transcripts for CD8_WT vs. non-CD8_WT"
	par(cex.main=0.8)
	sig_genes_de=sig_tn_de[,"id"]
	sig_gene_names_de=sig_tn_de[,"gene_name"]

	data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,labCol=short_names,col=rev(heat.colors(75)))


	dev.off()
	
	
