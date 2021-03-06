library(biomaRt)
library("dplyr")
library("intervals")

ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
attr = listAttributes(ensembl)

#exonData = getBM(c("ensembl_gene_id", "chromosome_name", "strand", "start_position", "end_position",
#               "ensembl_transcript_id", "transcript_start", "transcript_end",
#               "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", 
#               "gene_biotype", "external_gene_name"), mart=ensembl)
#exonData = filter(exonData, chromosome_name %in% c(1:19, "X", "Y"))

#exons = cbind(exonData$ensembl_exon_id, 
#              exonData$chromosome_name, exonData$exon_chrom_start, exonData$exon_chrom_end, 
#              exonData$ensembl_transcript_id, exonData$ensembl_gene_id)

exons = getBM(c("ensembl_exon_id", "chromosome_name", "strand", "exon_chrom_start", "exon_chrom_end",
                "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=ensembl)
#exons = filter(exons, chromosome_name %in% c(1:19, "X", "Y"))
exons = filter(exons, chromosome_name == 15)
exons = exons[order(exons$chromosome_name, exons$exon_chrom_start), ]
write.csv(exons, file = 'exons', sep = '\t', row.names=FALSE)

#exvect = rep(0, 104000000)
#by(exons[1:10,], 1:nrow(exons[1:10,]), function(item) {s = as.integer(item[4]); e = as.integer(item[5]);
                                         exvect[s:e] = exvect[s:e] + 1})
##apply(exons[1:10,], 1, function(item) {s = as.integer(item[4]); e = as.integer(item[5]);
##                                      exvect[s:e] = exvect[s:e] + 1})
#sum(exvect > 0)

#genes = getBM(c("ensembl_gene_id", "chromosome_name", "strand", "start_position", "end_position"), mart=ensembl)
              

snp_cast = read.table("CAST_15.bed")[,1:2]
snp_129S1 = read.table("129S1_15.bed")[,1:2]
colnames(snp_cast) = c("chromosome_name", "pos")
colnames(snp_129S1) = c("chromosome_name", "pos")
snp_cast$spec = "CAST"
snp_129S1$spec = "129S1"
snp = rbind(snp_cast, snp_129S1)
snp = snp[order(snp$chromosome_name, snp$pos), ]

snp_merged_vec = unique(snp$pos)
##snp_merged_count = sapply(snp_merged_vec, function(i) {sum(snp$pos == i)})
##snp_pos_count = cbind(snp_merged_vec, snp_merged_count)
#nonintersect = rep(0, 104000000)
#for (i in snp_cast$pos) {nonintersect[i] = nonintersect[i]+1}
#for (i in snp_129S1$pos) {nonintersect[i] = nonintersect[i]-1}
#nonintersect_snp = which(nonintersect != 0)

delta = 10000
x = c(0:(104000000%/%delta))*delta

y_cast = sapply(x, function(i) {sum(i <= snp_cast$pos & snp_cast$pos < i+delta)})
y_129S1 = sapply(x, function(i) {sum(i <= snp_129S1$pos & snp_129S1$pos < i+delta)})
plot (x, y_cast, type = 'h', ylim = c(-700,700), col = 'blue', 
      ylab = '# snp', xlab = 'position', main = paste('15, delta =',delta))
points (x, -y_129S1, type = 'h', col = 'red')

y_xor = sapply(x, function(i) {sum(i <= nonintersect_snp  & nonintersect_snp < i+delta)})
plot(x, y_xor, type = 'h')

y_unit = sapply(x, function(i) {sum(i <= snp_merged_vec & snp_merged_vec < i+delta)})
plot(x, y_unit, type = 'h', col = 'magenta',
     ylab = '# snp', xlab = 'position', main = paste('15, delta =',delta))


find_exon = function(pos, exondata) {
  res = exondata[exondata$exon_chrom_start <= pos &
                 exondata$exon_chrom_end >=pos,]
  if (nrow(res) == 0) {return (0)}
  else {print(pos); print(res);
        return (1)}
}
snp_region = snp_merged_vec[snp_merged_vec >= 7811054 & snp_merged_vec <= 8708742]
snp_region_exons = sapply(snp_region, function(i) {find_exon(i, exons)})
WOW = snp_region[which(snp_region_exons!=0)]
write.csv(cbind(rep(15, length(WOW)), WOW-1, WOW), file = 'exon_region_-1_snv', 
           row.names=FALSE)
           
chr15_snp_exons = sapply(snp_merged_vec[1:10400000], function(i) {find_exon(i, exons)})
snp_merged_vec[which(chr15_snp_exons != 0)]
snp_vs_exons = snp_merged_vec[which(chr15_snp_exons != 0)]  
write.table(cbind(rep(15, snp_vs_exons), snp_vs_exons-1, snp_vs_exons), file = 'exon_region_-1_snv', 
                                          row.names = FALSE, col.names = FALSE, sep = '\t')


y_exons = sapply(x, function(i) {sum(i <= WOW & WOW < i+delta)})
plot(x, y_exons, type = 'h', col = 'cyan',
     ylab = '# snp', xlab = 'position', main = paste('15, delta =',delta))


