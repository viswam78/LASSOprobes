#!/bin/bash

for i in *sorted;

do 

#../../../bowtie2-2.2.9/bowtie2 -x ../../k12CDS -U ${i} -S ${i}.sam -k 1 > ${i}.txt
#../../../samtools-1.3.1/samtools view -q 30 -bS ${i} > ${i}.bam
#../../../samtools-1.3.1/samtools sort ${i} > ${i}.sorted
#../../../bedtools2/bin/bedtools coverage -b ${i} -a ../../k12_chr_GN_targets.bed > targets${i}.out
#../../../bedtools2/bin/bedtools coverage -b ${i} -a ../../k12_chr_GN_nontargets.bed > nontargets${i}.out
#../../../bedtools2/bin/bedtools coverage -b ${i} -a ../../k12_chr_GN_complement.bed > complement${i}.out
#../../../bedtools2/bin/bedtools coverage -b ${i} -a ../../k12_chr_GN_complement.bed > complement${i}.out
../../../bedtools2/bin/genomeCoverageBed -ibam ${i} -g ../../k12_GN.bed -d > depthDGeneBG${i}.txt

done;
#../bowtie2-2.2.9/bowtie2-build k12/k12_chr.fna k12CDS
#../bowtie2-2.2.9/bowtie2 -x k12CDS -1 1_S1_L001_R1_001.fastq -2 1_S1_L001_R2_001.fastq -S HMVYVBCXX_1_CGATGT_1.sam -k 1 > orfAlignmentRate_k12_CDS.txt
#../../../bowtie2-2.2.9/bowtie2 -x ../../k12CDS -U HMVYVBCXX_1_CGATGT_1.fastq.gz -S HMVYVBCXX_1_CGATGT_1.sam -k 1 > HMVYVBCXX_1_CGATGT_1.txt
#../../../samtools-1.3.1/samtools view -q 30 -bS HMVYVBCXX_1_CGATGT_1.sam > HMVYVBCXX_1_CGATGT_1.bam
#../../../samtools-1.3.1/samtools sort HMVYVBCXX_1_CGATGT_1.bam > HMVYVBCXX_1_CGATGT_1.sorted
#../samtools-1.3.1/samtools mpileup -uf k12/k12_chr.fna HMVYVBCXX_1_CGATGT_1.sorted > HMVYVBCXX_1_CGATGT_1.sorted.mpileup
#../samtools-1.3.1/samtools index k12_CDS.sorted > k12_CDS.sorted.bai
#../bedtools2/bin/genomeCoverageBed -ibam k12_CDS.sorted -g k12/k12.ffn > orfHMVYVBCXX_1_CGATGT_12_CDS.bed
#../bedtools2/bin/bedtools multicov -bams k12_CDS.sorted -bed k12_chr_GN_targets.bed > orfCounts_k12_CDS_targets.txt
#../bedtools2/bin/bedtools multicov -bams k12_CDS.sorted -bed k12_Complement.bed > orfCountsComplement_k12_CDS.txt
#../bedtools2/bin/bedtools nuc -fi k12/k12.ffn -bed orfCounts_k12_CDS.txt > GC_k12_CDS.txt
#../bedtools2/bin/genomeCoverageBed -i k12_CDS.bed -g k12_chr_GN.bed -d > test.out
#../bedtools2/bin/bedtools coverage -b k12_CDS.sorted -a k12_chr_GN_targets.bed > test.out
#../bedtools2/bin/bedtools nuc -fi k12/k12.ffn -bed orfCountsComplement_k12_CDS.txt > GC_Complement_k12_CDS.txt

#../../../bedtools2/bin/bedtools coverage -b HMVYVBCXX_1_CGATGT_1.sorted -a ../../k12_chr_GN_targets.bed > targetsHMVYVBCXX_1_CGATGT_1.out
#../../../bedtools2/bin/bedtools coverage -b HMVYVBCXX_1_CGATGT_1.sorted -a ../../k12_chr_GN_nontargets.bed > nontargetsHMVYVBCXX_1_CGATGT_1.out
#../../../bedtools2/bin/bedtools coverage -b HMVYVBCXX_1_CGATGT_1.sorted -a ../../k12_chr_GN_complement.bed > complementHMVYVBCXX_1_CGATGT_1.out
#../bedtools2/bin/bedtools nuc -fi k12/k12_chr.fna -bed k12_chr_GN_targets.bed > targetsCoverageATGC.txt
#../bedtools2/bin/bedtools nuc -fi k12/k12_chr.fna -bed k12_chr_GN_nontargets.bed > nontargetsCoverageATGC.txt
#../bedtools2/bin/bedtools nuc -fi k12/k12_chr.fna -bed k12_chr_GN_complement.bed > complementCoverageATGC.txt
