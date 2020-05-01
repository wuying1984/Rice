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

### Raw data/resources:
#### *1: Genome scaffolds*: Bfra_R1V1.fa (reformat of consensus_MECAT_R1_90x.fasta From Dr. Menjun Hu)
#### *2: Botrytis cinerea EST sequence from NCBI (from Chris): EST_botryotinia_1.fasta
#### *3: Protein data: uniprot + all available Botrytis protein (uniprot_sprot.fasta 04182020; Botrytis protein from NCBI)
#### *4: Mycelia transcriptome: RNAseq data  

## 1. De Novo Repeat Identification
### 构建数据库
```
#!/bin/bash
#SBATCH -n 16
#SBATCH --mem=32G
cd /home/ywu/Xiaolab/Botrytis/Maker/RepeatModeler2
###Build database##################################
BuildDatabase -name Bfra ../Bfra_R1V1.fa
### run RepeatModeler###############################################################
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

## 2. map RNAseq data (from pure mycelia sample) to genome by hisat2 and assembly EST by trinity
### hisat2
```
~/Xiaolab/Botrytis/Maker/hisat$
hisat2-build Bfra_R1V1.fa Bfra_R1V1.fa 1>hisat2-build.log 2>&1
hisat2 --new-summary -p 6 -x Bfra_R1V1.fa -1 Bfra_1.fq  -2 Bfra_2.fq  2>hisat.log | samtools view -Sb - | samtools sort -o Bfra.bam -
```
### Trinity
upload *fastq.gz to web site https://galaxy.ncgas-trinity.indiana.edu/ (can not add --jaccard_clip)
OR (I Prefere to use the following)
```
#!/bin/bash
#SBATCH -n 32
#SBATCH --mem=64G
cd ~/Xiaolab/Botrytis/Maker/RNAseq
Trinity --seqType fq --left Bfra_1.fq --right Bfra_2.fq --jaccard_clip --max_memory 32G --CPU 24 --output trinity_out
```
### BRAKER training AUGUSTUS Gene prediction model
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


## 3. do maker
gffread my.gff3 -T -o my.gtf
#### genome BUSCO
```
run_BUSCO.py -i ../Bfra_R1V1.fa -l ~/program/BUSCO/ascomycota_odb9 -o Bfra_genome_asco  -m geno -c 6 -sp botrytis_cinerea >Bfra_genome_asco.out
```

### 1) First round
```
gff3_merge -d Bfra_R1V1.fa_master_datastore_index.log
fasta_merge -d Bfra_R1V1.fa_master_datastore_index.log
run_BUSCO.py -i Bfra_R1V1.fa.all.maker.proteins.fasta  -l ~/program/BUSCO/ascomycota_odb9 -m prot -c 4 -o protein.busco >protein.busco.out
```
####### C:96.8%[S:96.7%,D:0.1%],F:2.7%,M:0.5%,n:1315
####### INFO	1272 Complete BUSCOs (C)
####### INFO	1271 Complete and single-copy BUSCOs (S)
####### INFO	1 Complete and duplicated BUSCOs (D)
####### INFO	36 Fragmented BUSCOs (F)
####### INFO	7 Missing BUSCOs (M)
####### INFO	1315 Total BUSCO groups searched

### 2)Second round
# export 'confident' gene models from MAKER first round and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../Bfra_R1V1.fa_master_datastore_index.log
rename genome Bfra_R1V1_len50_aed0.25 genome.*
##### gather some stats and validate
`fathom Bfra_R1V1_len50_aed0.25.ann Bfra_R1V1_len50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom Bfra_R1V1_len50_aed0.25.ann Bfra_R1V1_len50_aed0.25.dna  -validate > validate.log 2>&1`
##### collect the training sequences and annotations, plus 1000 surrounding bp for training
`fathom Bfra_R1V1_len50_aed0.25.ann Bfra_R1V1_len50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1`
##### create the training parameters
```
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl Bfra_R1V1_R3 params > Bfra_R1V1_R3.hmm
```
```
gff3_merge -d Bfra_R1V1.fa_master_datastore_index.log
fasta_merge -d Bfra_R1V1.fa_master_datastore_index.log
run_BUSCO.py -i Bfra_R1V1.fa.all.maker.proteins.fasta  -l ~/program/BUSCO/ascomycota_odb9 -m prot -c 4 -o protein.busco.second >protein.busco.second.out
```

### 3)Third round
```
########GeneMark Training
########using protein######
perl /home/ywu/program_genome/gmes_linux_64/gmes_petap.pl  --seq ../Bfra_R1V1.fa.masked --EP --dbep ../../Evidence_files/Botrytis_protein/Bcin_B05.faa --fungus --verbose --cores=8 --soft_mask 2000 --min_contig 500 --max_intron 2000 1>gmes.e

########using only genome########
perl /home/ywu/program_genome/gmes_linux_64/gmes_petap.pl  --seq ../Bfra_R1V1.fa.masked --ES --fungus --verbose --cores=8 --soft_mask 2000 --min_contig 500 --max_intron 2000 1>gmes.e
```

```
gff3_merge -d Bfra_R1V1.fa_master_datastore_index.log
fasta_merge -d Bfra_R1V1.fa_master_datastore_index.log
run_BUSCO.py -i Bfra_R1V1.fa.all.maker.proteins.fasta  -l ~/program/BUSCO/ascomycota_odb9 -m prot -c 4 -o protein.busco.3rd >protein.busco.3rd.out
```

git status

git add *
git commit -m m 
git push

git pull
