#Adapter trimming
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -j 0 -o d9-5_trim.fq d9-5_FRRN210233347-1a_raw.fq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -j 0 -o dp16_trim.fq dp16_FRRN210233346-1a_raw.fq.gz

#Filter Rfam sequences 
##download Rfam sequences
grep "Gene; tRNA" Rfam-type.txt > tRNA.id
grep "Gene; rRNA" Rfam-type.txt > rRNA.id 
grep "snoRNA" Rfam-type.txt > snoRNA.id
grep "snRNA" Rfam-type.txt | grep -v "snoRNA"> snRNA.id
cat *id |cut -f 1 > need_download
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RFXXXXX.fa.gz
for i in  `cat tRNA.id |cut -f 1`;do cat fasta/${i}.fa >> Rfam.tRNA.fa; done
for i in  `cat rRNA.id |cut -f 1`;do cat fasta/${i}.fa >> Rfam.rRNA.fa; done
for i in  `cat snRNA.id |cut -f 1`;do cat fasta/${i}.fa >> Rfam.snRNA.fa; done
for i in  `cat snoRNA.id |cut -f 1`;do cat fasta/${i}.fa >> Rfam.snoRNA.fa; done
##make bowite index
bowtie-build Rfam.rRNA.fa Rfam.rRNA
bowtie-build Rfam.tRNA.fa Rfam.tRNA
bowtie-build Rfam.snoRNA.fa Rfam.snoRNA
bowtie-build Rfam.snRNA.fa Rfam.snRNA
##bowtie mapping and keep unmapped reads
sh filter_rfam.sh dp16_trim.fq
sh filter_rfam.sh dp9_5_trim.fq
cat filter_rfam.sh
bowtie -x ~/rfam/Rfam.rRNA --threads 15 ${1}.fq -S ${1}.sam --un ${1}.rRNA_.fq
bowtie -x ~/rfam/Rfam.tRNA --threads 15 ${1}.rRNA_.fq -S ${1}.sam --un ${1}.rRNA_tRNA_.fq
bowtie -x ~/rfam/Rfam.snoRNA --threads 15 ${1}.rRNA_tRNA_.fq -S ${1}.sam --un ${1}.rRNA_tRNA_snoRNA_.fq
bowtie -x ~/rfam/Rfam.snRNA --threads 15 ${1}.rRNA_tRNA_snoRNA_.fq -S ${1}.sam --un ${1}.rfam.filtered.fq

#Mapping to genome
sh get_specific_length.sh dp16
sh get_specific_length.sh d9-5
cat get_specific_length.sh:
seqkit seq -m 21 -M 22 ../1_rfam_filter/${1}.rfam.filtered.fq > ${1}_21_22.fq
seqkit seq -m 24 -M 24 ../1_rfam_filter/${1}.rfam.filtered.fq > ${1}_24.fq
sh mapping.sh dp16
sh mapping.sh d9-5
cat mapping.sh
bowtie -x ~/genome/tair10 --threads 15 -n 2 -a --best --strata ../1-2_trim_fastq/${1}_21_22.fq -S ${1}_21_22.sam >> genome_mapping_record 2>&1;
samtools view -bSF 4 ${1}_21_22.sam > ${1}_21_22.us.bam;
samtools sort -@ 10 -o ${1}_21_22.bam ${1}_21_22.us.bam;samtools index ${1}_21_22.bam;rm ${1}_21_22.us.bam
bowtie -x ~/genome/tair10 --threads 15 -n 2 -a --best --strata ../1-2_trim_fastq/${1}_24.fq -S ${1}_24.sam >> genome_mapping_record 2>&1;
samtools view -bSF 4 ${1}_24.sam > ${1}_24.us.bam;
samtools sort -@ 10 -o ${1}_24.bam ${1}_24.us.bam;samtools index ${1}_24.bam;rm ${1}_24.us.bam

#count the reads
featureCounts --largestOverlap -g ID -M -O -t gene -a TAIR10_GFF3_genes_transposons.gff -T 10 -o Novogene_TAIR10_read_count ../2.mapping/*.bam
featureCounts -M -O -T 10 -t gene -g Note -a gene.gff -o gene_stas ../2.mapping/*.bam
featureCounts -M -O -T 10 -t gene -g Alias -a te.gff -o te_stas ../2.mapping/*.bam
