<img src="./stacks_logo.png" width="30%" height="30%">

*Created by Narcis, July 2024, updated July 2024*

After finishing the A-tutorial, start with this one.

## Step 1: aligne the reads to the refrence genome using [bwa](https://bio-bwa.sourceforge.net)

### First, indext the genome

Example command:

```
bwa index -p Hpal /srv/ref_genome/out_JBAT_review3.FINAL_hap2.fa

```

### Second, align the samples to the refernce genome, convert SAM to BAM and sort the alignment




## Step 2: running [ref_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) program

Example command:


## Step 3: running [gstacks](https://catchenlab.life.illinois.edu/stacks/comp/gstacks.php)

Example command:



## Step 4: running [populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) program 

Example command:

```
populations -P dir -O dir -M popmap (filters) --fstats -k --sigma=150000 --bootstrap -N 100 --write-random-snp --structure --vcf --phylip --plink --genepop 
```









