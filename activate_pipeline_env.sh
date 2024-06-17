#!/bin/bash 
module load samtools/1.12  
module load java/jdk-17.0.2
module load bcftools/1.12 # this includes python version python3/3.9.2
source /g/data/pq84/python-3.9.2/bin/activate
module load bbmap/38.93
module load fastqc/0.11.7
#module load multiqc/1.9 
module load minimap2/2.28 # depen of RNA-Bloom
module load trinity/2.12.0
module load jellyfish/2.3.0 # trinity
module load bowtie2 # trinity
module load salmon/1.1.0 # trinity
export PATH=$PATH:/g/data/pq84/bin/RSEM-1.3.3 # trinity
export PATH=$PATH:/g/data/pq84/andrew/software/transabyss/:/g/data/pq84/andrew/software/blat/bin # transabyss, downstream of rna - bloom
module load abyss/2.2.3 # transabyss, downstream of rna - bloom
export PATH=$PATH:/g/data/pq84/bin/TransDecoder-TransDecoder-v5.7.1:/g/data/pq84/andrew/software/cdhit/ # downstream of rna - bloom
export PATH=$PATH:/g/data/pq84/bin/busco/bin:/g/data/pq84/bin/hmmer-3.4/bin:/g/data/pq84/bin/metaeuk/bin # BUSCO & its dependencies
module load R/4.3.1 # BUSCO
export R_LIBS=/home/588/jw1542/.local/lib/R_4.3:$R_LIBS #BUSCO
export PATH=$PATH:/g/data/pq84/bin/kallisto/build/bin 

export PATH=$PATH:/g/data/pq84/bin/corset-1.09-linux64 # maybe, but probably not neededpwd