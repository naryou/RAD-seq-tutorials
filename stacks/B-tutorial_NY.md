<img src="./stacks_logo.png" width="30%" height="30%">

*Created by Narcis, July 2024, updated August 2024*

After finishing the A-tutorial, start with this one.

## Step 1: align the reads to the refrence genome using [bwa](https://bio-bwa.sourceforge.net)

### First, indext the genome (create database)

Example command:

```
#!/bin/bash

module load Aligner/BWA/0.7.17
bwa index -p Hpal /srv/ref_genome/out_JBAT_review3.FINAL_hap2.fa

```

### Second, align the samples to the refernce genome, convert SAM to BAM and sort the alignment

```
#!/bin/bash

module load Aligner/BWA/0.7.17
module load Tools/samtools/1.17

src=/srv/kenlab/narcis/analysis/ddRAD/stacks
bwa_db=$src/bwa_index/Hpal
    
files=”SW_AG_Aa_1p
SW_BE_Lo_1t
SW_ZH_bu_1t
SW_ZH_We_1p
Fr_ELB_1p
Fr_ELB_4t
Fr_ELB_5t
Ge_MH5_1p
Ge_MH20_2p
Ge_MH21_1t
”

#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files
do 
    bwa mem -t 8 $bwa_db $src/samples/${sample}.fq.gz |
      samtools view -b |
      samtools sort --threads 4 > $src/bwa/aligned/${sample}.bam
done

#
```


## Step 2: running [ref_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) program (this step runs the whole pipeline, you can skip this step and run the next steps manually)

Example command:
```
ref_map.pl -T 8 --popmap ./popmaps/popmap -o ./ref_map/ --samples ./bwa/aligned
```

## Step 3: running [gstacks](https://catchenlab.life.illinois.edu/stacks/comp/gstacks.php)

Example command:

```
#!/bin/bash

module load Variants/Stacks/2.54

# Run gstacks to build loci from the aligned paired-end data.
#
gstacks -I $src/bwa/aligned/ -M $src/popmaps/popmap -O $src/stacks/ref_map -t 8
#
```


## Step 4: running [populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) program 

Example command:

```
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics and 
# smooth the statistics across the genome. Export several output files.
#
#!/bin/bash

module load Variants/Stacks/2.54

populations -P $src/stacks/ref_map -M $src/popmaps/popmap -r 0.65 --vcf --genepop --structure --fstats --smooth --hwe -t 8
#```









