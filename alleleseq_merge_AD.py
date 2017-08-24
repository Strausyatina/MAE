'''
Stephan Castel's code, modified:
-- add paired-end function
-- reduce memory burden (no sam uploading & chr divided)
-- in progress : merge counts  = STAR_qual - 2*(SNPmismatch) + 2*(SNPconsensus) - numOfAlignedOutOfExons

USAGE: ./alleleseq_merge_AD.py --pat_sam [path to file] --mat_sam [path to file] --o [path to folder] --paired [0|1]

'''

import sys;
import pandas;
import argparse;
import random;
import subprocess;
import time
start_time = time.time()

def load_sam_scores(sam, paired):
	print("LOADING SCORES SAM %s"%(sam))
	out_dict = {}
	in_stream = open(sam, "r")
	for line in in_stream:
        	if line[0] != "@": 
                	columns = line.strip().split("\t")
                	read_name = columns[0]
                	if int(columns[4]) != 255:
				continue
                    	## i.e. read is  uniquelly alligned 
                    	if (paired):
                        	if read_name not in out_dict:
                        		out_dict[read_name] = {}
                         	if int(columns[8]) > 0:
                             		out_dict[read_name]['A'] = get_score(line)
                         	else:
                             		out_dict[read_name]['B'] = get_score(line)
                   	else:
                        	out_dict[read_name] = get_score(line)
        in_stream.close()
	# NOTE: we suppose the data is correct, i.e. names are correct and unique, paired data consist of pairs , etc.
        return(out_dict)


def get_score(sam_line):
	q_score = 0
	sam_fields = sam_line.strip().split("\t")
	for field in sam_fields:
		if "AS:i:" in field: q_score = int(field.replace("AS:i:",""))
	return q_score


def main():
        #Arguments passed 
        parser = argparse.ArgumentParser()
        parser.add_argument("--pat_sam", required=True, help="Reads aligned to paternal genome")
        parser.add_argument("--mat_sam", required=True, help="Reads aligned to maternal genome")
        parser.add_argument("--o", required=True, help="Output file")
        parser.add_argument("--paired", default=0, help="Flag: If reads are paired-end")
        args = parser.parse_args()
        
        paired = int(args.paired)
        
        # First get reads that only map in one or the other
        pat_sam = load_sam_scores(args.pat_sam, paired)
        pat_reads = set(pat_sam.keys())
        mat_sam = load_sam_scores(args.mat_sam, paired)
        mat_reads = set(mat_sam.keys())

	mat_only_reads = mat_reads - pat_reads
        pat_only_reads = pat_reads - mat_reads

        # Find reads that are shared in both and chose the better alignment
        mat_count , pat_count = 0 , 0
        equal = 0
        mat_rand , pat_rand = 0 , 0
        mat_better_reads = set()
	pat_better_reads = set()
        
	print("MERGING READS");
        both_reads = mat_reads & pat_reads;
        for read in both_reads:
		if paired:
			mat_score , pat_score = mat_sam[read]['A']+mat_sam[read]['B'] , pat_sam[read]['A']+pat_sam[read]['B']
		else:  mat_score , pat_score = mat_sam[read], pat_sam[read]
                if mat_score > pat_score:
                        mat_count += 1
			mat_better_reads.add(read)
                elif pat_score > mat_score:
                        pat_count += 1;
                        pat_better_reads.add(read)
                else:
                        equal += 1;
                        x = random.randint(0, 1);
                        if x == 1:
				mat_better_reads.add(read)
                                mat_rand += 1
                        else:
				pat_better_reads.add(read)
                                pat_rand += 1

	# Get header
        header = subprocess.check_output("samtools view -SH "+args.mat_sam, shell=True)
        out_stream = open(args.o, "w")
        out_stream.write(header)
        out_stream.write("@RG     ID:pat\n")

	# Write reads 
	samfiles = [args.mat_sam, args.pat_sam]
	for samf in samfiles:
		if samf == args.mat_sam: samset = mat_better_reads.union(mat_only_reads) 
		else: samset = pat_better_reads.union(pat_only_reads)
        	in_stream = open(samf, "r")
        	for line in in_stream:
                	if line[0] != "@":
                        	if line.strip().split("\t")[0] in samset:
					 out_stream.write(line)
		in_stream.close()

        out_stream.close();        

        print("SUMMARY");
        print("%d reads in PAT only"%(len(pat_only_reads)));
        print("%d reads in MAT only"%(len(mat_only_reads)));
        print("%d reads where PATQ > MATQ"%(pat_count));
        print("%d reads where MATQ > PATQ"%(mat_count));
        print("%d reads where MATQ = PATQ"%(equal));
        print("%d MAT randomly selected"%(mat_rand));
        print("%d PAT randomly selected"%(pat_rand));
	print("----- %s seconds -----" % (time.time() - start_time))

if __name__ == "__main__":
        main();
