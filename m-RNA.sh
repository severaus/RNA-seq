#run trim_galore
mkdir trim_reads
trim_galore --quality 25 --length 70 --paired --fastqc -o ./trim_reads name_R1.fq.gz name_R2.fq.gz
#--quality 设定碱基质量的阈值，默认为20
#--length 设定输出reads的长度阈值，小于设定值会被抛弃
# --paired 对于双端测序结果，一对reads中，如果有一个被剔除，那么另一个也会被同样抛弃，而不管是否达到标准
# --fastqc 输入文件为fastqc文件
# -o 输出文件的路径

#build index
#genome index
hisat2-build -p 10 genome.fa genome #genome.fa代表参考基因组
# -p:设置线程
#annotation file
hisat2_extract_splice_sites.py genome.gtf>genome.splice_sites.bed #得到剪切位点的信息
hisat2_extract_exons.py genome.gtf>genome.exon.bed #得到外显子的信息
#以上2个python脚本是hisat2本身自带的
#hisat2——extract_snps.py snpcommon.txt>genoems.snp hiasat2也可以将SNP的信息加入到索引中
hisat2-build -p 10 --ss genome.splice_sites.bed --exon genome.exon.bed genome.fa genome 
#可以将基因组与注释文件一起建立索引，注释文件不是必须的，只是为了发现新的剪切位点

#mapping(hisat2+samtools+stringtie)
#hisat2比对
hisat2 -p 10 --dta -x /beegfs/project/uvm_mckay/WGBS-home/hisat2-cattle/ARS-UCD1.0.25 -1 name_R1.fq.gz -2 name_R2.fq.gz -S name.sam
# -p:线程数
# --dta:输出转录型的报告文件
# -x:基因组索引
# -S：输出sam文件
#samtools排序并将sam文件转化为bam文件
samtools sort -@ 8 -o name.bam name.sam 
# -@:设置线程数
# -o：输出文件


#!/bin/bash
#PBS -q batch


cd /storage2/zhangshengli/lwl/cattle/Cleandata_RNA/bzcj5
name="bzcj5"

export PATH=/apps/.local/software/bioinformatics/samtools-1.9/bin:$PATH

mkdir fastqc
fastqc --noextract -t 8 -f fastq -o ./fastqc ${name}_R1.fq.gz 
fastqc --noextract -t 8 -f fastq -o ./fastqc ${name}_R2.fq.gz 

hisat2 -p 10 --dta -x /storage2/zhangshengli/lwl/SRADATA/cattle/reference_genome/ARS-UCD1.2 -1 ${name}_R1.fq.gz -2 ${name}_R2.fq.gz -S ${name}.sam
samtools sort -@ 8 -o ${name}.bam ${name}.sam 
stringtie -p 10 -e -B -G /storage2/zhangshengli/lwl/SRADATA/cattle/reference_genome/ARS-UCD1.2.gtf -o ${name}.known.gtf -l ${name} ${name}.bam
