# CD8+ T cells
sk knockout vs. wildtype


## raw data processing


### check quality of reads

  fastqc *.fastq.gz 

  python3 -m multiqc . 

  mkdir fastqc
  mv *fastqc* fastqc
  
For mouse mm10 genome hisat2 pre-built index:

  mkdir -p $gbm/RNA_REF_FA
  cd $gbm/RNA_REF_FA
  wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
  tar -xzvf mm10_genome.tar.gz
