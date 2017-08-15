# input: 
#	PolyesterRun_uniformal.R [1]outputDir [2]faCAST [3]fa129S1 [4]paired|single [5]readlengths [6]mean_depth [7]%129S1 [8]error_rate 

print("Start Polyester simulating. Caslculating counts.")

args = commandArgs(trailingOnly=TRUE)

LOCALDIR = args[1] 
ENDS = args[4]
rlen = as.numeric(args[5])

setwd(LOCALDIR)

DEPTH = as.numeric(args[6])
Perc129S1 = as.numeric(args[7])
errr = as.numeric(args[8]) 

species = c('CAST', '129S1') 
faCAST = args[2]; fa129S1 = args[3] 
#"/home/amendelevich/Data/mice_genomes/pseudo/full15_naive/noErr_longestTr_uniform_fixedRPKM/CAST_15_cDNA_longestTr.fa"
#"/home/amendelevich/Data/mice_genomes/pseudo/full15_naive/noErr_longestTr_uniform_fixedRPKM/129S1_15_cDNA_longestTr.fa"

library("Biostrings") 

if (ENDS == 'paired') { pair = TRUE; print ("Simulating paired-end reads") }
else if (ENDS == 'single') { pair = FALSE; prinr ("Simulating single-end reads") } 
else {print ("ERROR in simulations, undefuned pairness of reads") 

###nTrs = 1103 
length129S1 = width(readDNAStringSet(fa129S1))
lengthCAST = width(readDNAStringSet(faCAST))
# counts = geneLen * coverage / (rlen * 2)
if (ENDS == 'paired') {
	Count129S1 = sapply(length129S1, function(x){as.integer(Perc129S1 * x * DEPTH / 100 / (rlen * 2)) +1})
	CountCAST = sapply(lengthCAST, function(x){as.integer((100-Perc129S1) * x * DEPTH/ 100 / (rlen * 2)) +1})
} else if (ENDS == 'single') {
        Count129S1 = sapply(length129S1, function(x){as.integer(Perc129S1 * x * DEPTH / 100 / (rlen)) +1})
        CountCAST = sapply(lengthCAST, function(x){as.integer((100-Perc129S1) * x * DEPTH/ 100 / (rlen)) +1})
} else {
	print("ERROR in simulations, undefuned pairness of reads")
}

library('polyester')

print("Simulating.")

simulate_experiment_countmat(fasta = fa129S1, readmat = t(t(Count129S1)),
    outdir = paste0(LOCALDIR,"/129S1"), 
    readlen = rlen, fragsd = 0, error_rate = errr, paired = pair)
simulate_experiment_countmat(fasta = faCAST, readmat = t(t(CountCAST)),
    outdir = paste0(LOCALDIR,"/CAST"), 
    readlen = rlen, fragsd = 0, error_rate = errr, paired = pair)
