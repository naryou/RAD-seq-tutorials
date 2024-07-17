<img src="./stacks_logo.png" width="30%" height="30%">

This is a tutorial to generate a VCF file using the [Stacks](https://catchenlab.life.illinois.edu/stacks/manual/) pipeline

Created by Narcis, July 2024, updated July 2024 



## Step 1: process_radtags

<img src="./process_radtags.png" width="70%" height="70%">



1)	Set up a directory structure on the server to hold your data in different steps

You do it as a “waterfall workspace”: working (clean (raw, samples), map (stacks))
Copy (don’t move!) the file containing raw Illumina sequences from “/home/workshop/Illumina_R1.fastq” into your raw directory

2)	Cleaning and de-multiplexing the sequences using process-radtags program manual [here](https://catchenlab.life.illinois.edu/stacks/manual/#clean)
You need to:
•	Specify a set of barcodes using `nano` or `vi`
Use these barcodes:

``` 
AACCC<tab>ind1
AATTT	ind2		
ACCAT	ind3	
CTCTT	ind4
GGCCT	ind5
```

•	Specify the restriction enzyme used to construct the library (SbfI)
•	Specify that process_radtags clean, discard, and rescue reads

Now, examine the results: the de-multiplexed sequences, and process_radtags.log file. What do you think about the content of the log file? Hint: use cat (for small files), head, more and tail.
Familiarize yourself with the format of the sequences (fastq format) and the quality (Phred scores). Here you find a good explanation: https://en.wikipedia.org/wiki/FASTQ_format 




## Step 2: denovo assembly
<img src="./denovo.png" width="75%" height="75%">

de novo assembly of RAD tags without a genome. Follow the manual [here](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php)

•	Run the Stacks’ denovo_map.pl program (*It’s important that you set –T (number of threads) to 1, so that you don’t occupy the whole server!)

o	You need to specify a map specifying which individuals belong to which population
Use `nano` or `vi` to make it as:

```
ind1<tab>1
ind2	1
ind3	1
ind4	2
ind5	2
```






## Step 3: population command



