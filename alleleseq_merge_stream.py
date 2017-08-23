##
## NOTE: 
## -- SAM-files should be already sorted by read names
##
## USAGE: alleleseq_merge_stream.py --pat_sam [path to file] --mat_sam [path to file] --o [path to file] --paired [0|1]
##

import sys
import pandas
import argparse
import random
import subprocess
import time
start_time = time.time()

def get_score(r_data, AS_pos, paired):
	q_score = 0
	if paired:
		reads_fields = [read.strip().split("\t") for read in r_data]
		q_score = int(reads_fields[0][AS_pos].replace("AS:i:",""))
	else:
		read_fields = r_data.strip().split("\t")
		q_score = int(read_fields[AS_pos].replace("AS:i:",""))
	return (q_score)

def get_AS(read, paired):
	if paired: read_fields = read[0].strip().split("\t")
	else: read_fields = read.strip().split("\t")
        for i in range(0, len(read_fields)):
        	if "AS:i:" in read_fields[i]: return (i)
	print ("No AS field")

def get_readname(data, paired):
	if paired: return (data[0].strip().split("\t")[0])
	else: return (data.strip().split("\t")[0])

def check_pairness(r_data):
        return (r_data[0].strip().split("\t")[0] == r_data[1].strip().split("\t")[0])

def check_readnames_eq(r_data1, r_data2, paired):
	return (get_readname(r_data1, paired) == get_readname(r_data2, paired))

def readline_to_write(r_data, paired):
        if paired: return (r_data[0] + r_data[1])
        return (r_data)

def main():

        # Parse arguments:

        parser = argparse.ArgumentParser()
        parser.add_argument("--pat_sam", required=True, help="Reads aligned to paternal genome")
        parser.add_argument("--mat_sam", required=True, help="Reads aligned to maternal genome")
        parser.add_argument("--o", required=True, help="Output file")
        parser.add_argument("--paired", default=0, help="Flag: If reads are paired-end")
        args = parser.parse_args()
        
        paired = int(args.paired)

        # Open output_sam; Get header:

        header = subprocess.check_output("samtools view -SH "+args.mat_sam, shell=True)
        out_stream = open(args.o, "w")
        out_stream.write(header)
        out_stream.write("@RG     ID:pat\n")

	# Open input_sam:        

        mat_count , mat_only , pat_count , pat_only = 0 , 0 , 0 , 0
        equal , mat_rand , pat_rand = 0 , 0 , 0

	source_m = open(args.mat_sam, 'r')
	source_p = open(args.pat_sam, 'r')

	bad_reads = set()

	def get_next_single(filehandler, prev_name):
		global bad_reads
		not_eq_to_last = True
		output = filehandler.readline()
		if output.strip().split('\t')[0] == prev_name: not_eq_to_last = False 
		while (bool(output) & (output.strip().split('\t')[0] == prev_name)):
			bad_reads.add(output.strip().split("\t")[0])
                        output = filehandler.readline()
        	return not_eq_to_last , output
	def get_next_pair(filehandler, prev_name):
		global bad_reads
                not_eq_to_last = True
		output = [filehandler.readline(), filehandler.readline()]
                if output[0].strip().split('\t')[0] == prev_name: not_eq_to_last = False
		while (bool(output[1]) & (output[0].strip().split('\t')[0] == prev_name) & (not check_pairness(output))):
			bad_reads.add(output[0].strip().split("\t")[0])
			output = [output[1], filehandler.readline()]
		if output[1]: return not_eq_to_last , output
		return not_eq_to_last , output[1]
	def get_data(filehandler, prev_name, paired):
		if paired: return(get_next_pair(filehandler, prev_name)) 
		else: return(get_next_single(filehandler, prev_name))

	# Skip header in each file:

	m_skip = int(subprocess.check_output("samtools view -SH "+args.mat_sam+" | wc -l", shell=True).strip())
	p_skip = int(subprocess.check_output("samtools view -SH "+args.pat_sam+" | wc -l", shell=True).strip())

	for i in range(m_skip):
		source_m.readline()
	for i in range(p_skip):
		source_p.readline()

	# Take first reads:
	m_read_name , p_read_name = '' , ''
	if paired:
		m_read , p_read = get_next_pair(source_m, m_read_name)[1] , get_next_pair(source_p, p_read_name)[1]
	else: m_read , p_read = get_next_single(source_m, m_read_name)[1] , get_next_single(source_p, p_read_name)[1]
	as_pos = get_AS(m_read, paired)

	# Merge till some EOF: 
	while (bool(m_read) & bool(p_read)):
		m_read_name , p_read_name = get_readname(m_read, paired) , get_readname(p_read, paired)
		m_next = get_data(source_m, m_read_name, paired) 
		p_next = get_data(source_p, p_read_name, paired)

		if m_next[0] & p_next[0] & check_readnames_eq(m_read, p_read, paired):
			# i.e. reads with equal names
			m_score = get_score(m_read, as_pos, paired)
                        p_score = get_score(p_read, as_pos, paired)
                        if m_score > p_score:
                                out_stream.write(readline_to_write(m_read, paired))
                                mat_count += 1
                        elif m_score < p_score:
                                out_stream.write(readline_to_write(p_read, paired))
                                pat_count +=1
                        else:
                                equal += 1
                                x = random.randint(0, 1)
                                if x == 1:
                                        mat_rand += 1
                                        out_stream.write(readline_to_write(m_read, paired))
                                else:
                                        pat_rand += 1
                                        out_stream.write(readline_to_write(p_read, paired))
			m_read , p_read = m_next[1] , p_next[1]
		elif m_next[0] & p_next[0] & (not check_readnames_eq(m_read, p_read, paired)):
			# i.e. unique reads
			if m_read_name > p_read_name : 
				pat_only += 1
				out_stream.write(readline_to_write(p_read, paired))
				p_read = p_next[1]
			else:
				mat_only += 1
				out_stream.write(readline_to_write(m_read, paired))
				m_read = m_next[1]
		else:
			# wrong reads
			if not m_next[0]: m_read = m_next[1]
			if not p_next[0]: p_read = p_next[1]
	
	print(pat_only, mat_only, equal, pat_count, mat_count)
	
	### Write the remain part of reads:
	if bool(m_read): 
		while bool(m_read):
			m_read_name = get_readname(m_read, paired)
			m_next = get_data(source_m, m_read_name, paired)
			if m_next[0]:
				mat_only += 1
                                out_stream.write(readline_to_write(m_read, paired))
                                m_read = m_next[1]
			else: m_read = m_next[1]
	elif bool(p_read):
		while bool(p_read):
			p_read_name = get_readname(p_read, paired)
			p_next = get_data(source_p, p_read_name, paired)
			if p_next[0]:
				pat_only += 1
                                out_stream.write(readline_to_write(p_read, paired))
                                p_read = p_next[1]
			else: p_read = p_next[1]
				

        out_stream.close(); source_m.close(); source_p.close()

        print("SUMMARY")
        print("%d reads in PAT only"%(pat_only))
        print("%d reads in MAT only"%(mat_only))
        print("%d reads where PATQ > MATQ"%(pat_count))
        print("%d reads where MATQ > PATQ"%(mat_count))
        print("%d reads where MATQ = PATQ"%(equal))
        print("%d MAT randomly selected"%(mat_rand))
        print("%d PAT randomly selected"%(pat_rand))
	print("----- %s seconds -----" % (time.time() - start_time))
	print("%d BAD read names: "%(len(bad_reads)) + " , ".join(sorted(list(bad_reads))))


if __name__ == "__main__":
        main();
