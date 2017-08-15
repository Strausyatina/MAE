#!/bin/bash

DTYPE="pseudo"
SIMDIR="~/MiceStats/Simulations/Pseudogenomes"
faCASTsim="/home/amendelevich/Data/mice_genomes/pseudo/full15_naive/noErr_longestTr_uniform_fixedRPKM/CAST_15_cDNA_longestTr.fa"
fa129S1sim="/home/amendelevich/Data/mice_genomes/pseudo/full15_naive/noErr_longestTr_uniform_fixedRPKM/129S1_15_cDNA_longestTr.fa"

for errr in 0 0.005; do
	for ps in single paired; do
		for depth in 10 25; do
			echo "./MiceStats.sh $DTYPE $SIMDIR $faCASTsim $fa129S1sim $ps 75 $depth 50 $errr"

done; done; done
