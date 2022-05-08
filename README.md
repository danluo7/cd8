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
    
    
##Index all bam files for IGV visualization

Make sure that the only files in the directory are the sam and bam/bai files, then:

can do: samtools index 1.bam

or:

    find *.bam -exec echo samtools index {} \; | sh


    

