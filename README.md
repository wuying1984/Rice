# Maker2_Botrytis_Fragariae_genome_annotion

## 1. De Novo Repeat Identification
### 构建数据库
```
/pub/software/RepeatModeler-2.0.1/BuildDatabase -name sesame test_data/genome.fasta
```

### 运行 RepeatModeler
`/pub/software/RepeatModeler-2.0.1/RepeatModeler -database sesame -pa 20 -LTRStruct`

### 运行 RepeatMasker
`/pub/software/RepeatMasker/RepeatMasker -pa 20 -qq -lib sesame-families.fa test_data/genome.fasta >repeatmasker.log 2>&1`

# 生成RepeatLandscape
/pub/software/RepeatMasker/util/calcDivergenceFromAlign.pl -s sesame.divsum test_data/genome.fasta.cat

perl /pub/software/RepeatMasker/util/createRepeatLandscape.pl -div sesame.divsum -g 18577337 > sesame.html

## 2. Full Repeat Annotation
First, we need to find the library name of Botrytis.
`cd ~/program/maker/exe/RepeatMasker2/RepeatMasker/util
./queryRepeatDatabase.pl -h`

First, we mask using a currated BovB/CR1 line library to overcome a previously-identified issue with the Repbase annotation. This probably won't be necessary in other species.

mkdir BovB_mask
RepeatMasker -pa 8 -e ncbi -lib bovb_cr1_species.lib -dir BovB_mask Boa_constrictor_SGA_7C_scaffolds.fa
Then the maksed FASTA from this search can be used as input for the next search, using the tetrapoda library from Repbase. I also normally rename the outputs after each round so they are more representative of what they contain.

mkdir Tetrapoda_mask
RepeatMasker -pa 8 -e ncbi -species tetrapoda -dir Tetrapoda_mask Boa_constrictor_SGA_7C.scaffolds.BovB.fa
And then I finished using two more rounds using a library of known and unknown snake repeats (including those from Boa). These rounds were split so that known elements would be preferentially annotated over unknown, to the degree possible.

mkdir 14Snake_known_mask 14Snake_unknown_mask
RepeatMasker -pa 8 -e ncbi -species 14Snakes_Known_TElib.fasta -dir 14Snake_known_mask Boa_constrictor_SGA_7C.scaffolds.Tetrapoda.fa
RepeatMasker -pa 8 -e ncbi -species 14Snakes_Known_TElib.fasta -dir 14Snake_known_mask Boa_constrictor_SGA_7C.scaffolds.14SnakeKnown.fa
Finally, results from each round must be analyzed together to produce the final repeat annotation

mkdir Full_mask
cp 14snake_unknown_mask/Boa_constrictor_SGA_7C.scaffolds.14SnakeUnknown.fa Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.fa
cp 14snake_unknown_mask/Boa_constrictor_SGA_7C.scaffolds.14SnakeUnknown.out Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out
gunzip BovB_mask/*.cat.gz Tetrapoda_mask/*.cat.gz 14snake_known_mask/*.cat.gz 14snake_unknown_mask/*.cat.gz
cat BovB_mask/*.cat Tetrapoda_mask/*.cat 14snake_known_mask/*.cat 14snake_unknown_mask/*.cat \
  > Final_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.cat
cd Full_mask
ProcessRepeats -species tetrapoda Final_mask/Boa_constrictor_SGA_7C.scaffolds.Full_mask.cat
Finally, in order to feed these repeats into MAKER properly, we must separate out the complex repeats (more info on this below).

# create GFF3
rmOutToGFF3.pl Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out > Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out.gff3
# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" Boa_constrictor_SGA_7C_scaffolds.full_mask.gff3 \
  > Boa_constrictor_SGA_7C_scaffolds.full_mask.complex.gff3
# reformat to work with MAKER
cat Boa_constrictor_SGA_7C_scaffolds.full_mask.complex.gff3 | \
  perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
  > Boa_constrictor_SGA_7C_scaffolds.full_mask.complex.reformat.gff3
Now we have the prerequesite data for running MAKER.
