<img src="./Structure.png" width="100%" height="100%">

This is the manual to create [Structure](https://web.stanford.edu/group/pritchardlab/structure.html) plots and visualize it using [CLUMPAK](https://tau.evolseq.net/clumpak/)

## Step 1: edit the STRUCTURE file 
The STRUCURE format file you have generated using the populations program in Stacks needs a small adjustment. 

* Remove the first line (line with a #)
* add a TAB (to have two TABs at the first line)

## Step 2: download the [Structure software](https://web.stanford.edu/group/pritchardlab/structure.html)

## Step 3: Run it as follows

* set K 1-6
* set burnin to 10000
* set MCMC to 20000

For the commandline version, edit the mainparams and extraparams files according to the input file. It's best to remove the marker names line and pop info column.

To get the number of loci, you need to count number of columns using this command : `cat file1 | awk 'BEGIN{FS=”\t”};{print NF}'`


## Step 4: Visualize the results in [CLUMPAK](https://tau.evolseq.net/clumpak/)

Select and zip the "f" files, and upload them in [CLUMPAK](https://tau.evolseq.net/clumpak/) 



