'''
Stephan Castel's code, modified:
-- add paired-end function

USAGE: 

NOTE:
1. SAM files should not contain multiple alignments (for this version) 
'''

import sys;
import pandas;
import argparse;
import random;
import subprocess;

def load_sam(sam, paired):
	print("LOADING SAM %s"%(sam))
	out_dict = {}
	in_stream = open(sam, "r")
	for line in in_stream:
        	if line[0] != "@":
                	columns = line.strip().split("\t")
                	read_name = columns[0]
                	if int(columns[4]) == 255:
                    		## i.e. read is  uniquelly alligned 
                    		if (paired):
                         		if read_name not in out_dict:
                             			out_dict[read_name] = {}
                         		if int(columns[8]) > 0:
                             			out_dict[read_name]['A'] = line
                         		else:
                             			out_dict[read_name]['B'] = line
                         	## it will work if all the paires are really pairs 
                         	## and we have no cases with "one-end pairs"
                   		else:
                        		out_dict[read_name] = line;
        in_stream.close()
        return(out_dict)


def get_star_alignment_score(sam_line, paired):
	if paired: 
		sam_line = sam_line['A']
        sam_fields = sam_line.strip().split("\t")
        for field in sam_fields:
                if "AS:i:" in field:
                        return(int(field.replace("AS:i:","")))
	return 0 


def readline_to_write(read_data, paired):
	if paired:
		return read_data['A']+read_data['B']
	return read_data


def main():
        #Arguments passed 
        parser = argparse.ArgumentParser()
        parser.add_argument("--pat_sam", required=True, help="Reads aligned to paternal genome")
        parser.add_argument("--mat_sam", required=True, help="Reads aligned to maternal genome")
        parser.add_argument("--o", required=True, help="Output file")
        parser.add_argument("--paired", default=0, help="Flag: If reads are paired-end")
        args = parser.parse_args()
        
        version = "0.1"
        print("")
        print("##################################################")
        print("              AlleleSeq Merge v%s"%(version))
        print("  Author: Stephane Castel (scastel@nygenome.org)")
        print("##################################################")
        print("")

        #1 load SAM files into memory
        paired = args.paired
        pat_sam = load_sam(args.pat_sam, paired)
        pat_reads = set(pat_sam.keys())
        mat_sam = load_sam(args.mat_sam, paired)
        mat_reads = set(mat_sam.keys())

        #1b get header
        header = subprocess.check_output("samtools view -SH "+args.mat_sam, shell=True)
        out_stream = open(args.o, "w")
        out_stream.write(header)
        out_stream.write("@RG     ID:pat\n")

        #2 first get reads that only map in one or the other
        mat_only_reads = mat_reads - pat_reads
        pat_only_reads = pat_reads - mat_reads
        for read in mat_only_reads:
		out_stream.write(readline_to_write(mat_sam[read], paired))
        for read in pat_only_reads:
                out_stream.write(readline_to_write(pat_sam[read], paired))

        #3 now find reads that are shared in both and chose the better alignment
        mat_count = 0;
        pat_count = 0;
        equal = 0;
        mat_rand = 0;
        pat_rand = 0;
        print("MERGING READS");
        both_reads = mat_reads & pat_reads;
        for read in both_reads:
                mat_score = get_star_alignment_score(mat_sam[read], paired)                
                pat_score = get_star_alignment_score(pat_sam[read], paired)
                if mat_score > pat_score:
                        mat_count += 1;
                        out_stream.write(readline_to_write(mat_sam[read], paired))
                
                elif pat_score > mat_score:
                        pat_count += 1;
                        out_stream.write(readline_to_write(pat_sam[read], paired))
                else:
                        equal += 1;
                        x = random.randint(0, 1);
                        if x == 1:
                                out_stream.write(readline_to_write(mat_sam[read], paired))
                                mat_rand += 1
                        else:
                                out_stream.write(readline_to_write(pat_sam[read], paired))
                                pat_rand += 1
        out_stream.close();        

        print("SUMMARY");
        print("%d reads in PAT only"%(len(pat_only_reads)));
        print("%d reads in MAT only"%(len(mat_only_reads)));
        print("%d reads where PATQ > MATQ"%(pat_count));
        print("%d reads where MATQ > PATQ"%(mat_count));
        print("%d reads where MATQ = PATQ"%(equal));
        print("%d MAT randomly selected"%(mat_rand));
        print("%d PAT randomly selected"%(pat_rand));


if __name__ == "__main__":
        main();
