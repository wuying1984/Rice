# MAKER2_PM_genome_annotation
## Using MAKER2 FOR ANNOTATION OF FOUR DICOT PM genomes
### Reference to  https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
###             & https://github.com/CompSynBioLab-KoreaUniv/FunGAP#step1

### Software & Data
### Software prerequisites:
###### RepeatModeler(2.0.1) and RepeatMasker (4.1.0) with all dependencies (I used NCBI BLAST) and RepBase (version used was 20181026).
###### MAKER version 2.31.9 (though any other version 2 releases should be okay).
###### Augustus version 3.3.3
###### BUSCO version 3.
###### SNAP https://github.com/KorfLab/SNAP
###### BEDtools version 2.24.0

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
RepeatMasker  -pa 8 -s -lib /home/ywu/Xiaolab/Botrytis/Maker/RepeatModeler2/Bfra-families.fa -nolow -xsmall Bfra_R1V1.fa >repeatmasker.log 2>&1
#####soft mask all reapeat
RepeatMasker  -pa 8 -s -lib /home/ywu/Xiaolab/Botrytis/Maker/RepeatModeler2/Bfra-families.fa -xsmall Bfra_R1V1.fa >repeatmasker.log 2>&1
```

# 生成RepeatLandscape
/pub/software/RepeatMasker/util/calcDivergenceFromAlign.pl -s sesame.divsum test_data/genome.fasta.cat

perl /pub/software/RepeatMasker/util/createRepeatLandscape.pl -div sesame.divsum -g 18577337 > sesame.html

# create GFF3
rmOutToGFF3.pl Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out > Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out.gff3


## 2. map RNAseq data to genome by hisat2 and assembly EST by trinity
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
#############use braker.pl in old version BRAKER######################################
perl /home/ywu/program_genome/BRAKER-2.1.2.old/scripts/braker.pl --species=BfraF --genome=../Bfra_R1V1.fa --bam=../hisat/Bfra.bam --prot_seq=/home/ywu/Xiaolab/Botrytis/Maker/Evidence_files/Botrytis_protein/Bcin_B05.faa --cores=8 --fungus --AUGUSTUS_CONFIG_PATH=/home/ywu/program/maker/exe/augustus/config/ --BAMTOOLS_PATH=/home/ywu/program_genome/software/bamtools/bin/ --SAMTOOLS_PATH=/home/ywu/program_genome/samtools --GENEMARK_PATH=/home/ywu/program_genome/gm_et_linux_64/gmes_petap --PYTHON3_PATH=/home/ywu/miniconda3/bin

braker.pl --species=Bfra --genome=../Bfra_R1V1.fa --bam=../hisat/Bfra.bam --prot_seq=../Bcin_B05.fasta  --cores=8 --AUGUSTUS_CONFIG_PATH=/home/ywu/program/maker/exe/augustus/config/ --BAMTOOLS_PATH=/home/ywu/program_genome/software/bamtools/bin/ --SAMTOOLS_PATH=/home/ywu/program_genome/samtools --GENEMARK_PATH=/home/ywu/program/maker/exe/gmes_petap

##############at gs8.genek.tv
##############use BRAKER-2.1.2
cd ~/workspace/Botrytis/BRAKER
braker.pl --species=BFRAR --genome=../Bfra_R1V1.fa.masked --softmasking   --bam=../Hisat/Bfra.bam --cores=8 --AUGUSTUS_CONFIG_PATH=/home/u2020/software/augustus-3.3.3/config/ --BAMTOOLS_PATH=/pub/software/bamtools/bin/ --SAMTOOLS_PATH=/pub/software/samtools --GENEMARK_PATH=/pub/software/gm_et_linux_64/gmes_petap

braker.pl --species=BFRAR --genome=../Bfra_R1V1.fa.masked --softmasking   --bam=../Hisat/Bfra.bam --fungus  --cores=8 --AUGUSTUS_CONFIG_PATH=/home/u2020/software/augustus-3.3.3/config/ --BAMTOOLS_PATH=/pub/software/bamtools/bin/ --SAMTOOLS_PATH=/pub/software/samtools --GENEMARK_PATH=/pub/software/gm_et_linux_64/gmes_petap

with protein
cd ~/workspace/Botrytis/BRAKER_protein
braker.pl --species=BFRAR --genome=../Bfra_R1V1.fa.masked --softmasking   --bam=../Hisat/Bfra.bam --prot_seq=../Bcin_B05.fasta --cores=8 --AUGUSTUS_CONFIG_PATH=/home/u2020/software/augustus-3.3.3/config/ --BAMTOOLS_PATH=/pub/software/bamtools/bin/ --SAMTOOLS_PATH=/pub/software/samtools --GENEMARK_PATH=/pub/software/gm_et_linux_64/gmes_petap
```

ll /pub/software/augustus-3.3.2/config/species/botrytis_cinerea/botrytis_cinerea
