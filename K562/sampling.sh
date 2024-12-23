#parallel bash sampling.sh ::: 0.04 0.06 0.08 
sample=$1
mkdir sample${sample}
zcat SRX15446041_S1_L001_R1_001.fastq.gz | seqkit sample -p ${sample} -s 10 -o sample${sample}/SRX15446041_S1_L001_R1_001.fastq.gz
zcat SRX15446041_S1_L001_R2_001.fastq.gz | seqkit sample -p ${sample} -s 10 -o sample${sample}/SRX15446041_S1_L001_R2_001.fastq.gz