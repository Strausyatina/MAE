# Pipeline, 2nd iteration

> Idea of changing method of obtaining data and using another simulator: Polyester

Given:
1. pseudo fasta for 2 lines (GRCm38 with distinguishing SNPs)
2. gene counts data from real experiments
---

1. Obtaining RNA-fasta for a region (15_7811054_8708742):
```
python3 mRNAfastaInterval_longestT.py --fasta genome.fa --gtf some.gtf --chrom N --s n1 --e n2 --ofasta mRNA_interval.fa --oGT longTrans_interval.tsv

```

2. Simulate RNA reads [Polyester]:
(https://github.com/leekgroup/polyester_code/blob/master/polyester_manuscript.Rmd)
* From given (by Sveta) gene counts obtain transcript counts (consider as longest transcript counts):
```
python3 countsIDsGtoT_allinone.py --GCMat SV1_AblY11_11_03_24_2015_R1_001_counts.txt --GenTrIDs longTrans_interval.tsv
```


---
And as before:

3. index references:
```
~/Tools/STAR --runMode genomeGenerate --genomeDir ~/Data/mice_genomes/pseudo/129S1_pseudo/ --genomeFastaFiles ~/Data/mice_genomes/pseudo/129S1_pseudo/129S1_pseudo.fa --sjdbGTFfile ~/Data/mice_genomes/pseudo/Mus_musculus.GRCm38.68.gtf
~/Tools/STAR --runMode genomeGenerate --genomeDir ~/Data/mice_genomes/pseudo/CAST_pseudo/ --genomeFastaFiles ~/Data/mice_genomes/pseudo/CAST_pseudo/CAST_pseudo.fa --sjdbGTFfile ~/Data/mice_genomes/pseudo/Mus_musculus.GRCm38.68.gtf
```

4. align simulated reads on the other reference:
```
~/Tools/STAR --genomeDir ~/Data/mice_genomes/pseudo/CAST_pseudo/ --readFilesIn ./129S1_pseudo_15_7811054_8708742-1.fastq ./129S1_pseudo_15_7811054_8708742-2.fastq --outFileNamePrefix ~/Data/mice_genomes/pseudo/129S1_pseudo_15_7811054_8708742.10000.

samtools view -Sb 129S1_pseudo_15_7811054_8708742.10000.Aligned.out.sam > 129S1_pseudo_15_7811054_8708742.10000.Aligned.out.bam
samtools sort 129S1_pseudo_15_7811054_8708742.10000.Aligned.out.bam 129S1_pseudo_15_7811054_8708742.10000
samtools index 129S1_pseudo_15_7811054_8708742.10000.bam

~/Tools/STAR --genomeDir ~/Data/mice_genomes/pseudo/129S1_pseudo/ --readFilesIn ./CAST_pseudo_15_7811054_8708742-1.fastq ./CAST_pseudo_15_7811054_8708742-2.fastq --outFileNamePrefix ~/Data/mice_genomes/pseudo/CAST_pseudo_15_7811054_8708742.10000.
  
samtools view -Sb CAST_pseudo_15_7811054_8708742.10000.Aligned.out.sam > CAST_pseudo_15_7811054_8708742.10000.Aligned.out.bam
samtools sort CAST_pseudo_15_7811054_8708742.10000.Aligned.out.bam CAST_pseudo_15_7811054_8708742.10000
samtools index CAST_pseudo_15_7811054_8708742.10000.bam
```

5. obtain bed file with snp positions for each pseudo_line:
```python
get_GTF_FASTA_sub.py obtainBEDfromVCF (vcf_name, bed_name)
```

6. mpileup:
```
samtools mpileup -I -B -q 0 -Q 0 -s -l 129S1_15.bed -f CAST_15pseudo/CAST_pseudo_15.fa 129S1_pseudo_15_7811054_8708742.10000.bam > 129S1_pileup
samtools mpileup -I -B -q 0 -Q 0 -s -l CAST_15.bed -f 129S1_15pseudo/129S1_pseudo_15.fa CAST_pseudo_15_7811054_8708742.10000.bam |head   
```
