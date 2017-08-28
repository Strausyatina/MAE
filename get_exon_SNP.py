# NOTE: need sort by coords and only exon lines from gtf
# if CHROM without 'chr':
# grep "	exon	" GRCm38.gtf | sort -k1n -k4,5n > /home/asya/Downloads/GRCm38_15_exons_sort.gtf
# NOTE: vcf should contain snp only

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp_vcf", required=True, help="VCF with only SNP")
    parser.add_argument("--gtf_exon", required=True, help="GTF with exons only")
    parser.add_argument("--o", required=True, help="Output file")
    args = parser.parse_args()

    out_stream = open(args.o, 'w')
    vcf = open(args.snp_vcf , 'r')
    gtf = open(args.gtf_exon , 'r')

    vcf_line = vcf.readline()
    while(vcf_line.startswith("##")):
        vcf_line = vcf.readline()
    vcf_header = vcf_line.replace('#','').strip().split()
    for i in range(len(vcf_header)):
        if vcf_header == "POS": pos_col = i 
	else if vcf_header == "CHROM": chrom_col = i

    def get_next_snp(vcf_handler):
        line = vcf_handler.readline()
        if line:
            line_cont = line.strip().split()
            return line_cont[chrom_col] , line_cont[pos_col] , line
        else: return ''
    def get_next_exon(gtf_handler):
        line_cont = gtf_handler.readline().strip().split()
        return line_cont[0] , line_cont[3] , line_cont[4]
    def snp_in_interval(SNP, I):
        snp_c , snp_p , _ = SNP
        I_c , I_s , I_e = I
        if I_c < snp_c or (I_c == snp_c and I_e < snp_p): return 1
        elif I_c > snp_c or (I_c == snp_c and I_s > snp_p): return -1
        elif I_s <= snp_p and snp_p <= snp_e: return 0

    snp_cur = get_next_snp(vcf)
    exon_cur = get_next_exon(gtf_handler)
    while snp_cur and exon_cur:
        sii = snp_in_interval(snp_cur, exon_cur)
        if sii == 0: out_stream.write(snp_cur[2])
        elif sii > 0: exon_cur = get_next_exon(gtf_handler)
        else: snp_cur = get_next_snp(vcf)  

    out_stream.close() ; vcf.close() ; gtf.close()


