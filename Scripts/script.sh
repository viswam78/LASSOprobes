#!/bin/bash

for i in trimmomatic*.gz.sam.bam.Sorted.bed;

do 

#java -jar ../Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 ${i} trimmomatic${i} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#../bowtie2-2.2-2.9/bowtie2 -x k12/k12NEW -p 4 --very-sensitive --no-unal -U ${i} -S ${i}.sam -k 1 > ${i}.txt
#../samtools-1.3.1/samtools view -q 30 -bS ${i} > ${i}.bam
#../samtools-1.3.1/samtools sort ${i} > ${i}.Sorted
#../bedtools2/bin/genomeCoverageBed -ibam ${i} -g k12/k12strands.bed -d > depthDGeneBG${i}.txt
#../bedtools2/bin/bamToBed -i ${i}  > ${i}.bed
../bedtools2/bin/coverageBed -b ${i} -a k12/targetk12chrhundred.bed -wa > targetshundred${i}.out
#../bedtools2/bin/coverageBed -b ${i} -a k12/nontargetk12chr.bed -wa > nontargets${i}.out
#../bedtools2/bin/coverageBed -b ${i} -a k12/complementk12chr.bed > complement${i}.out
../bedtools2/bin/coverageBed -b ${i} -a k12/targetk12chrhundred.bed  -mean -wa > targetshundredd${i}.out
#../bedtools2/bin/coverageBed -b ${i} -a k12/nontargetk12chr.bed  -mean -wa > nontargetsd${i}.out
#../bedtools2/bin/coverageBed -b ${i} -a k12/complementk12chr.bed  -mean > complementd${i}.out

done;

#../bedtools2/bin/bedtools nuc -fi k12/k12_chr.fna -bed k12/targetk12chr.bed > targetsCoverageATGC.txt
#../bedtools2/bin/bedtools nuc -fi k12/k12_chr.fna -bed k12/nontargetk12chr.bed > nontargetsCoverageATGC.txt
#../bedtools2/bin/bedtools nuc -fi k12/k12_chr.fna -bed k12/complementk12chr.bed > complementCoverageATGC.txt
