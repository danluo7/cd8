# CD8+ T cells
sk knockout vs. wildtype


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

take a look at the files to see the platform info etc:

        gzip -cd Wild-CDT3_S3_L001_R2_001.fastq | head

output: 
@A00257:731:HV577DRXY:1:2101:1289:1000 2:N:0:TATGCCTTAC+TAATGTGTCT
TCATGTTCAGACCTCACACTTGAAAAAGAAAATTTTCTACCCCCATGGTCC


### At this point decide on naming conventions. Order the samples logically

hisat2 -p 8 --rg-id=HV577DRXY.1 --rg SM:1 --rg LB:1_TATGCCTTAC+TAATGTGTCT --rg PL:ILLUMINA --rg PU:HV577DRXY.1.TATGCCTTAC+TAATGTGTCT -x $cd8/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_1_S1_L001_R1_001.fastq.gz -2 $gbm_data/6931_1_S1_L001_R2_001.fastq.gz -S $gbm_data/alignments/1_rep1.sam


