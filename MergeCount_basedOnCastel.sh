#!/bin/bash
# /home/amendelevich/Scripts/MergeCount_basedOnCastel.sh $NNAME $ALDIR $MERGEDIR $REF129S1 $REFCAST $GTF129S1 $GTFCAST $PS $VCF $MTYPE

NNAME=$1
ALDIR=$2
MERGEDIR=$3
REF129S1=$4
REFCAST=$5
GTF129S1=$6
GTFCAST=$7
PS=$8
VCF=$9
MTYPE=$10

if [[ $MTYPE = "simple" ]]; then
	if [[ $PS = "paired" ]]; then
		for t in pair left; do
			python ~/Tools/MAE/alleleseq_merge_UPD.py --paired 1 --pat_sam "$ALDIR/$NNAME/$NNAME."$t".onCAST.Aligned.out.sam" --mat_sam "$ALDIR/$NNAME/$NNAME."$t".on129S1.Aligned.out.sam" --o "$ALDIR/$NNAME/$NNAME."$t".merged.sam" > "$ALDIR/$NNAME/$NNAME."$t".merged.log"
            		/home/tools/samtools-1.2/samtools view -bS "$ALDIR/$NNAME/$NNAME."$t".merged.sam" | /home/tools/samtools-1.2/samtools sort -@ 10 -T $BASE -o "$ALDIR/$NNAME/$NNAME."$t".merged.bam" -
            		/home/tools/samtools-1.2/samtools index "$ALDIR/$NNAME/$NNAME."$t".merged.bam"
            		/home/tools/samtools-1.2/samtools idxstats "$ALDIR/$NNAME/$NNAME."$t".merged.bam" > "$ALDIR/$NNAME/$NNAME."$t".merged.idxstats"
		done
	else if [[ $PS = "single" ]]; then
                        python ~/Tools/MAE/alleleseq_merge_UPD.py --paired 1 --pat_sam "$ALDIR/$NNAME/$NNAME.single.onCAST.Aligned.out.sam" --mat_sam "$ALDIR/$NNAME/$NNAME.single.on129S1.Aligned.out.sam" --o "$ALDIR/$NNAME/$NNAME.single.merged.sam" > "$ALDIR/$NNAME/$NNAME.single.merged.log"
                        /home/tools/samtools-1.2/samtools view -bS "$ALDIR/$NNAME/$NNAME.single.merged.sam" | /home/tools/samtools-1.2/samtools sort -@ 10 -T $BASE -o "$ALDIR/$NNAME/$NNAME.single.merged.bam" -
                        /home/tools/samtools-1.2/samtools index "$ALDIR/$NNAME/$NNAME.single.merged.bam"
                        /home/tools/samtools-1.2/samtools idxstats "$ALDIR/$NNAME/$NNAME.single.merged.bam" > "$ALDIR/$NNAME/$NNAME.single.merged.idxstats"
	else echo "ERROR: bad ends type"
	fi
else if [[ MTYPE = "advanced" ]]; then
	#
	### the same with prefix _AD_ ?
	#
else echo "ERROR: bad merge type"
fi

##for item in merged 129S1 CAST; do
##/home/tools/samtools-1.2/samtools stats -c 1,100,1    .bam    | grep ^COV | cut -f 2- > $LOCALDIR"/statcounts_merged.singleend"
##            /home/tools/samtools-1.2/samtools stats -c 1,100,1 $BASE"_on129S1.singleend.Aligned.out.bam" | grep ^COV | cut -f 2- > $LOCALDIR"/statcounts_129S1.singleend"
##            /home/tools/samtools-1.2/samtools stats -c 1,100,1 $BASE"_onCAST.singleend.Aligned.out.bam" | grep ^COV | cut -f 2- > $LOCALDIR"/statcounts_CAST.singleend"

##python ~/Tools/MAE/allelecounter_samV.py --vsamtools /home/tools/samtools-1.2/samtools --vcf $VCF --sample F1 --bam $BASE".merged.singleend.bam" --ref ~/Data/mice_genomes/pseudo/129S1_15pseudo/129S1_pseudo_15.fa --min_cov 0 --min_baseq 0 --min_mapq 10 --o $LOCALDIR/129S1-ref_statcounts.singleend
### TO DO: с видоизменённым vcf, где вписаны доп.мутации.

