#!/bin/bash
#PBS -d .
#PBS -l walltime=10:00:00,mem=20gb

DTYPE="pseudo"
SIMDIR="~/MiceStats/Simulations/Pseudogenomes"
faCASTsim="/home/amendelevich/Data/mice_genomes/pseudo/full15_naive/noErr_longestTr_uniform_fixedRPKM/CAST_15_cDNA_longestTr.fa"
fa129S1sim="/home/amendelevich/Data/mice_genomes/pseudo/full15_naive/noErr_longestTr_uniform_fixedRPKM/129S1_15_cDNA_longestTr.fa"

dCASTref="/home/amendelevich/Data/mice_genomes/pseudo/CAST_15pseudo/"
d129S1ref="/home/amendelevich/Data/mice_genomes/pseudo/129S1_15pseudo/"
gtfCAST="/home/amendelevich/Data/mice_genomes/pseudo/GRCm38_15.gtf"
gtf129S1="/home/amendelevich/Data/mice_genomes/pseudo/GRCm38_15.gtf"
ALDIR="~/MiceStats/Alignments/Pseudogenomes"
MTYPE="simple"
MDIR="~/MiceStats/Merge/Pseudogenomes"


#       MiceStats.sh [1]typeofdata [2]SIMDir [3]faCAST [4]fa129S1 [5]paired|single [6]readlengths [7]mean_depth [8]%129S1 [9]error_rate 
#                    [10]referenceCAST [11]reference129S1 [12]ALDir [13]gtfCAST [14]gtf129S1
#                    [15]mergeType [16]mergeDir

for errr in 0 0.005; do
	for ps in single paired; do
		for depth in 10 25; do
			echo "./MiceStats.sh $DTYPE $SIMDIR $faCASTsim $fa129S1sim $ps 75 $depth 50 $errr $dCASTref $d129S1ref $ALDIR $gtfCAST $gtf129S1 $MTYPE $MDIR"
			~/Scripts/MiceStats.sh $DTYPE $SIMDIR $faCASTsim $fa129S1sim $ps 75 $depth 50 $errr $dCASTref $d129S1ref $ALDIR $gtfCAST $gtf129S1 $MTYPE $MDIR
done; done; done
