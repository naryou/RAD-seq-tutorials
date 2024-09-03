First step after getting your data from the sequencing center is to look at the quality if the data. We will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for that aim.

```
fastqc /home/course/raw_data/*.gz -o /home/course/YOURNAME/fastqc --extract --threads 2
```

Now inspect the generated plots. How does the quality look like?


