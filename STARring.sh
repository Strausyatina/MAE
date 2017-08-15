#!/bin/bash
# input:
#       STARring.sh $NNAME $SIMDIR $ALDIR $REF129S1 $REFCAST $GTF129S1 $GTFCAST $PS

NNAME=$1
SIMDIR=$2
ALDIR=$3
REF129S1=$4
REFCAST=$5
GTF129S1=$6
GTFCAST=$7
PS=$8

if [[ $PS = "paired" ]]; then
	sed -i 's/>read/>129S1read/g'  $SIMDIR/$NNAME/129S1/sample_01_*.fasta
	sed -i 's/>read/>CASTread/g'  $SIMDIR/$NNAME/CAST/sample_01_*.fasta
	cat $SIMDIR/$NNAME/129S1/sample_01_1.fasta $SIMDIR/$NNAME/CAST/sample_01_1.fasta > $SIMDIR/$NNAME/simreads_1.fasta
	cat $SIMDIR/$NNAME/129S1/sample_01_2.fasta $SIMDIR/$NNAME/CAST/sample_01_2.fasta > $SIMDIR/$NNAME/simreads_2.fasta

	# PAIRED:
	/home/amendelevich/Tools/STAR --genomeDir $REF129S1 --readFilesIn $SIMDIR/$NNAME/simreads_1.fasta $SIMDIR/$NNAME/simreads_2.fasta --sjdbGTFfile $GTF129S1 --outFileNamePrefix $ALDIR/$NNAME/$NNAME".pair.on129S1." --runThreadN 10 --outFilterMultimapNmax 1 --outSAMattrRGline ID:mat  --quantMode TranscriptomeSAM  --outSAMattributes All
	/home/amendelevich/Tools/STAR --genomeDir $REFCAST --readFilesIn $SIMDIR/$NNAME/simreads_1.fasta $SIMDIR/$NNAME/simreads_2.fasta --sjdbGTFfile $GTFCAST --outFileNamePrefix $ALDIR/$NNAME/$NNAME".pair.onCAST." --runThreadN 10 --outFilterMultimapNmax 1 --outSAMattrRGline ID:mat  --quantMode TranscriptomeSAM  --outSAMattributes All

	# LEFT END:
	~/Tools/STAR --genomeDir $REF129S1 --readFilesIn $SIMDIR/$NNAME/simreads_1.fasta --sjdbGTFfile $GTF129S1 --outFileNamePrefix $ALDIR/$NNAME/$NNAME".left.on129S1." --runThreadN 10 --outFilterMultimapNmax 1 --outSAMattrRGline ID:mat  --quantMode TranscriptomeSAM  --outSAMattributes All
	~/Tools/STAR --genomeDir $REFCAST --readFilesIn $SIMDIR/$NNAME/simreads_1.fasta --sjdbGTFfile $GTFCAST --outFileNamePrefix $ALDIR/$NNAME/$NNAME".left.onCAST." --runThreadN 10 --outFilterMultimapNmax 1 --outSAMattrRGline ID:mat  --quantMode TranscriptomeSAM  --outSAMattributes All

	# BAM from SAM :
        for item in 129S1 CAST; do
		for t in pair left; do 
        		/home/tools/samtools-1.2/samtools view -bS "$ALDIR/$NNAME/$NNAME."$t".on"$item".Aligned.out.sam" | /home/tools/samtools-1.2/samtools sort -@ 10 -T "$ALDIR/$NNAME/$NNAME."$t".on"$item -o "$ALDIR/$NNAME/$NNAME."$t".on"$item".Aligned.out.bam" -
        		/home/tools/samtools-1.2/samtools index "$ALDIR/$NNAME/$NNAME."$t".on"$item".Aligned.out.bam"
	done; done;

else if [[ $PS = "single" ]]; then
	sed -i 's/>read/>129S1read/g'  $SIMDIR/$NNAME/129S1/sample_01.fasta
        sed -i 's/>read/>CASTread/g'  $SIMDIR/$NNAME/CAST/sample_01.fasta
        cat $SIMDIR/$NNAME/129S1/sample_01.fasta $SIMDIR/$NNAME/CAST/sample_01.fasta > $SIMDIR/$NNAME/simreads.fasta

	# SINGLE :
	/home/amendelevich/Tools/STAR --genomeDir $REF129S1 --readFilesIn $SIMDIR/$NNAME/simreads.fasta --sjdbGTFfile $GTF129S1 --outFileNamePrefix $ALDIR/$NNAME/$NNAME".single.on129S1." --runThreadN 10 --outFilterMultimapNmax 1 --outSAMattrRGline ID:mat  --quantMode TranscriptomeSAM  --outSAMattributes All
        /home/amendelevich/Tools/STAR --genomeDir $REFCAST --readFilesIn $SIMDIR/$NNAME/simreads.fasta --sjdbGTFfile $GTFCAST --outFileNamePrefix $ALDIR/$NNAME/$NNAME".single.onCAST." --runThreadN 10 --outFilterMultimapNmax 1 --outSAMattrRGline ID:mat  --quantMode TranscriptomeSAM  --outSAMattributes All

	# BAM FROM SAM :
	for item in 129S1 CAST; do
        	/home/tools/samtools-1.2/samtools view -bS "$ALDIR/$NNAME/$NNAME.single.on"$item".Aligned.out.sam" | /home/tools/samtools-1.2/samtools sort -@ 10 -T "$ALDIR/$NNAME/$NNAME.single.on"$item -o "$ALDIR/$NNAME/$NNAME.single.on"$item".Aligned.out.bam" -
                /home/tools/samtools-1.2/samtools index "$ALDIR/$NNAME/$NNAME.single.on"$item".Aligned.out.bam"
        done;

else echo "ERROR: type should be paired or single" 
fi


