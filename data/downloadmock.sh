#!/bin/sh
#SBATCH -p rome
#SBATCH -t 03:00:00

# Downloads the Mock 5 and Mock 12 reads, metadata and reference sequences from Mockrobiota.

#Change into the mock community numbers to download 
mock="5 12"

# change directory to where to download the files
cd /gpfs/work2/0/gs/danique/mockrobiota/

# used data inventory: https://github.com/caporaso-lab/mockrobiota/blob/master/inventory.tsv
# Loop through mock community numbers 
for num in $mock; do

 echo $num

 # create new directory for mock file if it does not exist yet
 if ! [ -d "mock$num" ]; then
   mkdir "mock$num"
   mkdir "mock$num/seq"
   mkdir "mock$num/seq/output"
   mkdir "mock$num/info"
   mkdir "mock$num/fastq_demux"
 fi
 
 # raw sequences
 cd "mock$num/seq"

 wget "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-$num/mock-forward-read.fastq.gz"
 wget "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-$num/mock-reverse-read.fastq.gz"
 wget "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-$num/mock-index-read.fastq.gz"
 
 mv mock-forward-read.fastq.gz forward.fastq.gz
 mv mock-reverse-read.fastq.gz reverse.fastq.gz
 mv mock-index-read.fastq.gz barcodes.fastq.gz
 
 # sample info
 cd "../info"
  
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/README.md"
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/dataset-metadata.tsv"
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/sample-metadata.tsv"
 
 # silva 99% otu
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/silva/123/99-otus/database-identifiers.tsv"
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/silva/123/99-otus/expected-taxonomy.tsv"
 
 mv database-identifiers.tsv silva-database-identifiers.tsv
 mv expected-taxonomy.tsv silva-expected-taxonomy.tsv
 
 # source data 
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/source/expected-sequences.fasta"
 wget "https://github.com/caporaso-lab/mockrobiota/raw/master/data/mock-$num/source/taxonomy.tsv"
 
 cd /gpfs/work2/0/gs/danique/mockrobiota/

done