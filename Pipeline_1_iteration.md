# Pipeline, 1st iteration

> This pipeline is bad: we do not need CDS, only exons; wgsim does something strange around snps. 

1. from FASTA files obtain "RNA" FASTA files with the longest CDS "transcripts" for genes on the interval:

```python
get_GTF_FASTA_sub.py obtainCHRfromFASTA (fasta, chrFasta, chromosome) 
```
obtain fasta with only 15 chromosome
```python      
get_GTF_FASTA_sub.py getGTFsubGeneCDSonSegm (gtf, output_name, output_name_filtered, I) 
```
obtain subGTF for interval
```
sort filtered_file
```
```python
get_GTF_FASTA_sub.py obtainSubFasta (fasta_file, I_file, output_file, header):
```      
obtain "RNA" fasta

2. simulate RNA reads:
```
~/Tools/wgsim/wgsim -1 70 -2 70 -N 10000  129S1_pseudo_15__7811054_8708742.fa 129S1_pseudo_15_7811054_8708742-1.fastq 129S1_pseudo_15_7811054_8708742-2.fastq > out-129S1-10000
~/Tools/wgsim/wgsim -1 70 -2 70 -N 10000  CAST_pseudo_15__7811054_8708742.fa CAST_pseudo_15_7811054_8708742-1.fastq CAST_pseudo_15_7811054_8708742-2.fastq > out-129S1-10000
```

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
