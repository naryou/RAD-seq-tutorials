First step after getting your data from the sequencing center is to look at the **integrity** and **quality** of the data. For integrity we use [MD5Sum](https://en.wikipedia.org/wiki/Md5sum#:~:text=md5sum%20is%20used%20to%20verify,error%20or%20non%2Dmalicious%20meddling.) check and for quality We will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Example command for integrity check:

```
nohup md5sum *.gz > md5sum.txt &  
nohup md5sum -c md5sum.txt md5sum_NY.txt > md5check.txt &

```


Example command for quality check:

```
fastqc /home/course/raw_data/*.gz -o /home/course/YOURNAME/fastqc --extract --threads 2
```

Now inspect the generated plots. How does the quality look like?


