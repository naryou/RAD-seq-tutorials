First step after getting your data from the sequencing center is to look at the **integrity** and **quality** of the data. For integrity we use [MD5Sum](https://en.wikipedia.org/wiki/Md5sum#:~:text=md5sum%20is%20used%20to%20verify,error%20or%20non%2Dmalicious%20meddling.) check and for quality We will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

When we recieve data from genomic centers or other sources, there is a file called `md5sum.txt` accompanying it. We want to generate the same file and compare it with the source file to make sure that the size of the downloaded file is the same as the source file. 

Example command for integrity check:

```
nohup md5sum *.gz > md5sum_YOUR_INITIALS.txt &  
nohup md5sum -c md5sum.txt md5sum_YOUR_INITIALS.txt > md5check.txt &

```


Example command for quality check:

```
nohup fastqc /home/course/raw_data/*.gz -o /home/course/YOURNAME/fastqc --extract --threads 2 &
```

Now inspect the generated plots. How does the quality look like?


