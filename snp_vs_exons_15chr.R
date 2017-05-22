## for one chromosome only (here -- 15)

library(biomaRt)
library("dplyr")

ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")

# obtain exon data for 15 chr:

exons = getBM(c("ensembl_exon_id", "chromosome_name", "strand", "exon_chrom_start", "exon_chrom_end",
                "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=ensembl)
exons = filter(exons, chromosome_name == 15)
exons = exons[order(exons$chromosome_name, exons$exon_chrom_start), ]

# obtain snp data for 15 chr:

snp_cast = read.table("CAST_15.bed")[,1:2]; colnames(snp_cast) = c("chromosome_name", "pos")
snp_129S1 = read.table("129S1_15.bed")[,1:2]; colnames(snp_129S1) = c("chromosome_name", "pos")
snp_cast$spec = "CAST"; snp_129S1$spec = "129S1"
snp = rbind(snp_cast, snp_129S1)
snp = snp[order(snp$chromosome_name, snp$pos), ]

snp_merged_vec = unique(snp$pos)

# obtain "this snp leave in exon" vector (it's ... slow =)) :

find_exon = function(pos, exondata) {
  res = exondata[exondata$exon_chrom_start <= pos &
                 exondata$exon_chrom_end >=pos,]
  if (nrow(res) == 0) {return (0)}
  else {return (1)}
}

# desired (testing) part:
#snp_region = snp_merged_vec[snp_merged_vec >= 7811054 & snp_merged_vec <= 8708742]
#snp_region_exons = sapply(snp_region, function(i) {find_exon(i, exons)})
#WOW = snp_region[which(snp_region_exons!=0)]
#write.csv(cbind(rep(15, length(WOW)), WOW-1, WOW), file = 'exon_region_-1_snv', 
#           row.names=FALSE)
## (-1) because of mpileup...

chr15_snp_exons = sapply(snp_merged_vec, function(i) {find_exon(i, exons)})
snp_vs_exons = snp_merged_vec[which(chr15_snp_exons != 0)]

# bed file for mpileup:

write.table(cbind(rep(15, snp_vs_exons), snp_vs_exons-1, snp_vs_exons), file = 'exon_region_-1_snv', 
           row.names = FALSE, col.names = FALSE, sep = '\t')
## (!) BAM should be sorted

##--------------------------------------------------------------------------------------------------------------------------
## distribution picts (each spec. and united):

delta = 10000
x = c(0:(104000000%/%delta))*delta

y_cast = sapply(x, function(i) {sum(i <= snp_cast$pos & snp_cast$pos < i+delta)})
y_129S1 = sapply(x, function(i) {sum(i <= snp_129S1$pos & snp_129S1$pos < i+delta)})
plot (x, y_cast, type = 'h', ylim = c(-700,700), col = 'blue', 
      ylab = '# snp', xlab = 'position', main = paste('15, delta =',delta))
points (x, -y_129S1, type = 'h', col = 'red')

y_unit = sapply(x, function(i) {sum(i <= snp_merged_vec & snp_merged_vec < i+delta)})
plot(x, y_unit, type = 'h', col = 'magenta',
     ylab = '# snp', xlab = 'position', main = paste('15, delta =',delta))
