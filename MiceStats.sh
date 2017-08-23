#!/bin/bash
# input:
#	MiceStats.sh [1]typeofdata [2]SIMDir [3]faCAST [4]fa129S1 [5]paired|single [6]readlengths [7]mean_depth [8]%129S1 [9]error_rate 
#		     [10]referenceCAST [11]reference129S1 [12]ALDir [13]gtfCAST [14]gtf129S1
#		     [15]mergeType [16]mergeDir

SIMDIR=$2
TYPEDATA=$1
fa129S1sim=$4
faCASTsim=$3
PS=$5
RLEN=$6
MDEPTH=$7
PROP=$8
ERRR=$9

NAME=$TYPEDATA"_129S1_"$PROP"_CAST_md"$MDEPTH"_"$PS"_l"$RLEN"_err"$ERRR

echo "Begin 10 simulations for $NAME :"
for N in 0 1 2 3 4 5 6 7 8 9; do
    SIMNAME=$NAME"_"$N 
    NSIMDIR=$SIMDIR/$SIMNAME
    mkdir -p $NSIMDIR
    echo "Rscript /home/amendelevich/Scripts/PolyesterRun_uniformal.R $NSIMDIR $faCASTsim $fa129S1sim $PS $RLEN $MDEPTH $PROP $ERRR"
    Rscript /home/amendelevich/Scripts/PolyesterRun_uniformal.R $NSIMDIR $faCASTsim $fa129S1sim $PS $RLEN $MDEPTH $PROP $ERRR
done

#REFCAST=$10
#REF129S1=$11
#ALDIR=$12
#GTFCAST=$13
#GTF129S1=$14

#echo "Begin aligning by 10 simulated fastas by STAR :"
#for N in 0 1 2 3 4 5 6 7 8 9; do
#    NNAME=$NAME"_"$N
#    mkdir -p $ALDIR/$NNAME
#    /home/amendelevich/Scripts/STARring.sh $NNAME $SIMDIR $ALDIR $REF129S1 $REFCAST $GTF129S1 $GTFCAST $PS 
#done

#MTYPE=$15
#MERGEDIR=$16

#echo "Begin allele_merge and allele_counting :"
#for N in 0 1 2 3 4 5 6 7 8 9; do
#    NNAME=$NAME"_"$N
#    mkdir $MERGEDIR/$NNAME
#    /home/amendelevich/Scripts/MergeCount_basedOnCastel.sh $NNAME $ALDIR $MERGEDIR $REF129S1 $REFCAST $GTF129S1 $GTFCAST $PS $VCF $MTYPE 
#done
