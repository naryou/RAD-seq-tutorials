#!/bin/bash
src1=/home/ubuntu/course/
src2=$src1/raw_data

$src1/executables/hts_SeqScreener --singleend-input $src2/R1_001.fastq.gz \
-L R1_001.phix.json -f R1_001.phiX 
