# Maker2_Botrytis_Fragariae_genome_annotion

## 1. De Novo Repeat Identification
### 构建数据库
```
#!/bin/bash
#PBS -l ncpus=16
#PBS -l mem=32G
cd /home/ywu/Xiaolab/Botrytis/Maker/RepeatModeler2
BuildDatabase -name Bfra ../Bfra_R1V1.fa
```

### 运行 RepeatModeler
```
RepeatModeler -database Bfra -pa 16 -LTRStruct 1>repeatmodeler2.o 2>repeatmodeler2.e
```

### 运行 RepeatMasker
```
RepeatMasker -pa 16 -s -lib /home/ywu/Xiaolab/Botrytis/Maker/RepeatModeler2/Bfra-families.fa Bfra_R1V1.fa >repeatmasker.log 2>&1
#####soft mask + no masking of low complexity
RepeatMasker  -pa 8 -s -lib /home/ywu/Xiaolab/Botrytis/Maker/RepeatModeler2/Bfra-families.fa -nolow -small Bfra_R1V1.fa >repeatmasker.log 2>&1
```

# 生成RepeatLandscape
/pub/software/RepeatMasker/util/calcDivergenceFromAlign.pl -s sesame.divsum test_data/genome.fasta.cat

perl /pub/software/RepeatMasker/util/createRepeatLandscape.pl -div sesame.divsum -g 18577337 > sesame.html

# create GFF3
rmOutToGFF3.pl Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out > Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out.gff3


## 2. map RNAseq data to genome by hisat2 and trinity
```
~/Xiaolab/Botrytis/Maker/hisat$
hisat2-build Bfra_R1V1.fa Bfra_R1V1.fa 1>hisat2-build.log 2>&1
hisat2 --new-summary -p 6 -x Bfra_R1V1.fa -1 Bfra_1.fq  -2 Bfra_2.fq  2>hisat.log | samtools view -Sb - | samtools sort -o Bfra.bam -
```
```
upload *fastq.gz to web site https://galaxy.ncgas-trinity.indiana.edu/ (can not add --jaccard_clip)
OR (I Prefere to use the following)
#!/bin/bash
#PBS -l ncpus=16
#PBS -l mem=32G
cd ~/Xiaolab/Botrytis/Maker/RNAseq
Trinity --seqType fq --left Bfra_1.fq --right Bfra_2.fq --jaccard_clip --max_memory 32G --CPU 24 --output trinity_out
```

## BRAKER training 
```
braker.pl --species=Bfra --genome=../Bfra_R1V1.fa --bam=../hisat/Bfra.bam --cores-8 --AUGUSTUS_CONFIG_PATH=/home/ywu/program/maker/exe/augustus/config/ --BAMTOOLS_PATH=/home/ywu/program_genome/software/bamtools/bin/ --SAMTOOLS_PATH=/home/ywu/program_genome/samtools --GENEMARK_PATH=/home/ywu/program/maker/exe/gmes_petap
```

