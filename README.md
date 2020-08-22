# MAKER2_Rice_genome_annotation
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
#### *1: Genome scaffolds*: From Dr. Zhou
#### *2: rice EST sequence from rice paper (japonica_3transcripts.fasta, rufipongon_3transcripts.fasta, nivara_3transcripts.fast); EST from PlantGDB rice_EST_from_PlantGDB.fasta;     used file /scratch/Xiao_Group/Rice/Evidence/trinity/3species_plantGDB.fasta 
#### *3: Protein data: uniprot + all available rice protein (uniprot_sprot.fasta 04182020; http://ibi.zju.edu.cn/ricerelativesgd/download.php) used file /scratch/Xiao_Group/Rice/Evidence/Protein/uni_addOryza.fasta_20200503

## 1. De Novo Repeat Identification
### 构建数据库及运行RepeatMasker
```
#!/bin/bash
#SBATCH -n 16
#SBATCH --mem=32G
cd /scratch/Xiao_Group/Rice/RepeatModeler
echo ZJ1
echo ZJ1_Genome_HERA.fasta
mkdir ZJ1
cd ZJ1
cp $i ${i/_Genome_HERA.fasta/}
BuildDatabase -name ZJ1 ../ZJ1_Genome_HERA.fasta
RepeatModeler -database ZJ1 -pa 16 -LTRStruct
#####RepeatMasker from now#############
cd /scratch/Xiao_Group/Rice/
mkdir ZJ1_RepeatMasker
cd ZJ1_RepeatMasker
mv /scratch/Xiao_Group/Rice/RepeatModeler/ZJ1_Genome_HERA.fasta ./
RepeatMasker -pa 16 -s -lib /scratch/Xiao_Group/Rice/RepeatModeler/ZJ1/ZJ1-families.fa -xsmall ZJ1_Genome_HERA.fasta >ZJ1_repeatmasker.log 2>&1
```

## 2. do maker
split genome
```
cat  ZJ5_Genome_HERA.masked.fa | g ">" | sed 's#>##' >list.txt
for i in `cat list.txt`; do cat ZJ5_Genome_HERA.masked.fa | grep -A1 "$i$" >part.$i.fa; done
cat part.Contig*.fa >Chr13.fa
rename part.Chr Chr part.Chr*
rename .fa .masked.fa *.fa
```
```
for i in Chr*.fa ; do echo $i; cat maker_opts.ctl | sed "s#genome_file#$i#"  >maker_opts.ctl_${i/.masked.fa/};done
```

creat file "do_maker.sh"
```
#!/bin/bash
#SBATCH -n 16
#SBATCH --mem=32G
cd /scratch/Xiao_Group/Rice/ZJ5
maker opt maker_bopts.ctl maker_exe.ctl --fix_nucleotides -base result 2>oe/error
```
```
for i in maker_opts.ctl_*;do echo $i;cat do_maker.sh| sed "s#optctl#$i#;s#result#${i/maker_opts.ctl/ZJ5}#;s#error#${i/maker_opts.ctl/ZJ5_R1}.e#" >ZJ5_${i/maker_opts.ctl_/}_maker.sh;done
```
for i in {1..4} ; do echo $i;echo ZJ5_Chr${i}_maker.sh;done ##USE 16 cpu and 32G mem
for i in {5..13} ; do echo $i;echo ZJ5_Chr${i}_maker.sh;done ##USE 8 cpu and 12G mem
gffread my.gff3 -T -o my.gtf

#### genome BUSCO
ZJ5_BUSCO.sh
```
#!/bin/bash
#SBATCH -n 16
#SBATCH --mem=32G
cd /scratch/Xiao_Group/Rice/BUSCO_genome
echo ZJ5
echo ZJ5_Genome_HERA.fasta
mkdir ZJ5
cd ZJ5
run_BUSCO.py -i /scratch/Xiao_Group/Rice/RepeatModeler/ZJ5_Genome_HERA.fasta -l ~/program/BUSCO/liliopsida_odb10 -o ZJ5_genome_lili -m geno -c 16 -sp maize >ZJ5_genome_lili.out
run_BUSCO.py -i /scratch/Xiao_Group/Rice/RepeatModeler/ZJ5_Genome_HERA.fasta -l ~/program/BUSCO/liliopsida_odb10 -o ZJ5_genome_lili_long -m geno -c 16 -sp maize --long >ZJ5_genome_lili_long.out
```
####### C:90.6%[S:89.5%,D:1.1%],F:4.1%,M:5.3%,n:3278

####### C:95.8%[S:94.7%,D:1.1%],F:2.1%,M:2.1%,n:3278
####### INFO	3140 Complete BUSCOs (C)
####### INFO	3103 Complete and single-copy BUSCOs (S)
####### INFO	37 Complete and duplicated BUSCOs (D)
####### INFO	69 Fragmented BUSCOs (F)
####### INFO	69 Missing BUSCOs (M)


### 1) First round
for i in {01..12} ; do echo $i; cat $i/$i.fasta|grep -v ">" | fold -w3000000 >${i}_split_3M;done
wc -l *split_3M
for i in {1..15}; do echo $i; echo ">01_$i" >01_${i}.fasta; cat 01_split_3M | line $i >>01_${i}.fasta;done


15 01_split_3M
       13 02_split_3M
       14 03_split_3M
       13 04_split_3M
       11 05_split_3M
       11 06_split_3M
       11 07_split_3M
       10 08_split_3M
        9 09_split_3M
        9 10_split_3M
       12 11_split_3M
       10 12_split_3M

```
for i in ZJ5_Chr*.maker.output; do echo $i; cd $i; gff3_merge -d *_index.log; fasta_merge -d *_index.log;cd ../;done
for i in ZJ5_Chr*.maker.output; do echo $i; cd $i; cp *.gff ../ROUND1_result; cp *.fasta ../ROUND1_result; cd ../;done
cat ZJ5_Chr*.all.maker.proteins.fasta >ZJ5_R1.all.maker.proteins.fasta
run_BUSCO.py -i ZJ5_R1.all.maker.proteins.fasta  -l ~/program/BUSCO/liliopsida_odb10 -m prot -c 4 -o ZJ5_R1_PROTEIN >ZJ5_R1_PROTEIN.busco.out

for i in {1..15};do echo $i;echo ">$i" >ZJ4_01_${i}.fasta;cat split_3M.fasta| line $i >>ZJ4_01_${i}.fasta;done
for i in {1..15} ; do echo $i; cat maker_opts_ZJ4_01.ctl| sed "s#01\_1#01\_${i}#" >maker_opts_ZJ4_01_${i}.ctl; done
for i in {1..15} ; do echo $i; cat ZJ4_01_01_ROUND1.sh| sed "s#01\_01#01\_${i}#g" >01_${i}.sh; done
for i in ZJ4_01_0*.maker.output; do echo $i; cd $i; gff3_merge -d *_index.log; fasta_merge -d *_index.log;cd ../;done
for i in ZJ4_01_0*.maker.output; do echo $i; cd $i; cp *.gff ../ROUND1_result_2; cp *.all.maker.proteins.fasta ../ROUND1_result_2; cd ../;done

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
