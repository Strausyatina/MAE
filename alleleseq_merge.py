'''
by Stephan Castel
'''

import sys;
import pandas;
import argparse;
import random;
import subprocess;

def main():
	#Arguments passed 
	parser = argparse.ArgumentParser()
	parser.add_argument("--pat_sam", help="Reads aligned to paternal genome")
	parser.add_argument("--mat_sam", help="Reads aligned to maternal genome")
	parser.add_argument("--o", help="Output file")
	args = parser.parse_args()
	
	version = "0.1";
	print("");
	print("##################################################")
	print("              AlleleSeq Merge v%s"%(version));
	print("  Author: Stephane Castel (scastel@nygenome.org)")
	print("##################################################");
	print("");
	
	#1 load SAM files into memory
	pat_sam = load_sam(args.pat_sam);
	pat_reads = set(pat_sam.keys());
	mat_sam = load_sam(args.mat_sam);
	mat_reads = set(mat_sam.keys());

	
	#1b get header
	header = subprocess.check_output("samtools view -SH "+args.pat_sam, shell=True).replace("_paternal","");
	
	out_stream = open(args.o, "w");
	out_stream.write(header);
	out_stream.write("@RG     ID:mat\n");
	
	#2 first get reads that only map in one or the other
	pat_only_reads = pat_reads - mat_reads;
	for read in pat_only_reads:
		out_stream.write(pat_sam[read])
	
	mat_only_reads = mat_reads - pat_reads;
	for read in mat_only_reads:
		out_stream.write(mat_sam[read])
	
	#3 now find reads that are shared in both and chose the better alignment
	mat_count = 0;
	pat_count = 0;
	equal = 0;
	mat_rand = 0;
	pat_rand = 0;	
	print("MERGING READS");
	both_reads = mat_reads & pat_reads;
	for read in both_reads:
		pat_score = get_star_alignment_score(pat_sam[read])
		mat_score = get_star_alignment_score(mat_sam[read])
		if pat_score > mat_score:
			pat_count += 1;
			out_stream.write(pat_sam[read])
		elif pat_score < mat_score:
			mat_count += 1;
			out_stream.write(mat_sam[read])
		else:
			equal += 1;
			x = random.randint(0, 1);
			if x == 1:
				out_stream.write(pat_sam[read])
				pat_rand += 1
			else:
				out_stream.write(mat_sam[read])
				mat_rand += 1
	
	out_stream.close();
	
	print("SUMMARY");
	print("%d reads in PAT only"%(len(pat_only_reads)));
	print("%d reads in MAT only"%(len(mat_only_reads)));
	print("%d reads where PATQ > MATQ"%(pat_count));
	print("%d reads where MATQ > PATQ"%(mat_count));
	print("%d reads where MATQ = PATQ"%(equal));
	print("%d MAT randomly selected"%(mat_rand));
	print("%d PAT randomly selected"%(pat_rand));	

def load_sam(sam):
	print("LOADING SAM %s"%(sam));
	out_dict = {};
	in_stream = open(sam, "r");
	for line in in_stream:
		if line[0:1] != "@":
			#READNAME						   MAPQ CIGAR
			#ERR188053.20303311	163	1	13785	255	75M	=	13897	177	CAGGAGGCTGCCATTTGTCCTGCCCACCTTCTTAGAAGCGAGACGGAGCAGACCCATCTGCTACTGCCCTTTCTA	CCCFFFFFHHHHHJJJJJJJJJIJIIJJJJJJJJJJJIIJJJIIJJJJJJHJIJJJHHHHHFFFFFFEDEEDDDD	PG:Z:MarkDuplicates	RG:Z:id	NH:i:1	HI:i:1	nM:i:0	AS:i:138
			#ERR188053.20303311	83	1	13897	255	65M10S	=	13785	-177	CTAGTCTCAATTTAAGAAGATCCCCATGGCCACAGGGCCCCTGCCTGGGGGCTTGTCACCTCCCCACCTTCTTCC	DDDDDDDDEEDEDDDDDDDBBDDEFFDFFHHAHCJJIJJJJJIJJJJJIIIJJJJIHEJIIJHHFAHFFFFFBCB	PG:Z:MarkDuplicates	RG:Z:id	NH:i:1	HI:i:1	nM:i:0	AS:i:138
			
			sam_line = line.replace("_maternal","").replace("_paternal","");
			columns = sam_line.split("\t");
			read_name = columns[0];
			if int(columns[4]) == 255:
				out_dict[read_name] = sam_line;
			#	if int(columns[8]) > 0:
			#		out_dict[read_name+'.A'] = sam_line;
			#	else:
			#		out_dict[read_name+'.B'] = sam_line;
	in_stream.close()
	
	return(out_dict);	

def get_star_alignment_score(sam_line):
	sam_fields = sam_line.replace("\n","").split("\t");
	for field in sam_fields:
		if "AS:i:" in field:
			return(int(field.replace("AS:i:","")));

if __name__ == "__main__":
	main();
