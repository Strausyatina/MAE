##
## NOTE: 
## -- SAM-files should be already sorted by read names
## -- Read names must be non-empty
## USAGE: alleleseq_merge_stream.py --pat_sam [path to file] --mat_sam [path to file] --o [path to file] --paired [0|1]
##

import sys
import pandas
import argparse
import random
import subprocess
import time
start_time = time.time()

# special samtools ordering:
def isdigit(c):
    return c.isdigit()
def uporotiycmp(firstrb, secstrb):
    firstr = firstrb + '\x00' 
    secstr = secstrb + '\x00'
    firind , secind = 0 , 0
    firlen = len(firstr)
    seclen = len(secstr)
    while (firind < firlen) and (secind < seclen):
        if (isdigit(firstr[firind]) and isdigit(secstr[secind])):
            while (firstr[firind] == '0'):
                firind += 1
            while (secstr[secind] == '0'):
                secind += 1
            while isdigit(firstr[firind]) and isdigit(secstr[secind]) and (firstr[firind] == secstr[secind]):
                firind += 1
                secind += 1
            if (isdigit(firstr[firind]) and isdigit(secstr[secind])):
                i = 0
                while (isdigit(firstr[firind+i]) and isdigit(secstr[secind+i])):
                    i += 1
                if isdigit(firstr[firind+i]): return 1
                elif isdigit(secstr[secind+i]): return -1
                else: return ord(firstr[firind]) - ord(secstr[secind])
            elif (isdigit(firstr[firind])): return 1
            elif (isdigit(secstr[secind])): return -1
            elif (firind < secind): return 1
            elif (firind > secind): return -1
        else:
            if (firstr[firind] != secstr[secind]):
               return ord(firstr[firind]) - ord(secstr[secind])
            firind += 1
            secind += 1
    if firind < firlen: return 1
    elif secind < seclen: return -1
    else: return 0

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

def get_name(line):
    if not line: return ''
    return line.strip().split('\t')[0]
def get_readname(data, paired):
    if paired: return get_name(data[0])
    else: return get_name(data)

def readline_to_write(r_data, paired):
    if paired: return (r_data[0] + r_data[1])
    return (r_data)

def asc_order(name1, name2):
    return (uporotiycmp(name1, name2) <= 0)


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
    bad_reads = set()

    source_m = open(args.mat_sam, 'r')
    source_p = open(args.pat_sam, 'r')

    def get_next_block(fhandler):
        first_line = fhandler.readline()
        if first_line == '': return [] , ''
        first_name = get_name(first_line)
        output = []
        cur_line = first_line 
        while get_name(cur_line) == first_name:
            output.append(cur_line)
            cur_pos = fhandler.tell()
            cur_line = fhandler.readline()
        fhandler.seek(cur_pos)
        return output , first_name
    def get_next_single(fhandler):
        cur_block = get_next_block(fhandler)
        while len(cur_block[0]) > 1:
            bad_reads.add(cur_block[1])
            cur_block = get_next_block(fhandler)
        if not cur_block[0]: return ''
        return cur_block[0][0]
    def get_next_pair(fhandler):
        cur_block = get_next_block(fhandler)
        while (len(cur_block[0]) == 1) or (len(cur_block[0]) > 2):
            bad_reads.add(cur_block[1])
            cur_block = get_next_block(fhandler)
        if not cur_block[0]: return []
        return cur_block[0]
        
    # Skip header in each file:

    m_skip = int(subprocess.check_output("samtools view -SH "+args.mat_sam+" | wc -l", shell=True).strip())
    p_skip = int(subprocess.check_output("samtools view -SH "+args.pat_sam+" | wc -l", shell=True).strip())

    for i in range(m_skip):
        source_m.readline()
    for i in range(p_skip):
        source_p.readline()

    def update_m():
        if paired: return get_next_pair(source_m)
        else: return get_next_single(source_m) 
    def update_p():
        if paired: return get_next_pair(source_p)
        else: return get_next_single(source_p) 

    # Take first reads:

    m_read , p_read = update_m() , update_p() 
    as_pos = get_AS(m_read, paired)

    # Merge till some EOF:
 
    while (bool(m_read) and bool(p_read)):
        m_read_name , p_read_name = get_readname(m_read, paired) , get_readname(p_read, paired)
        if m_read_name == p_read_name:
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
            m_read , p_read = update_m() , update_p()
        else:
            if asc_order(p_read_name, m_read_name):
                pat_only += 1
                out_stream.write(readline_to_write(p_read, paired))
                p_read = update_p()
            else:
                mat_only += 1
                out_stream.write(readline_to_write(m_read, paired))
                m_read = update_m()


    print(pat_only, mat_only, equal, pat_count, mat_count)
   
 
    ### Write the remain part of reads:
    if bool(m_read): 
        while bool(m_read):
            mat_only += 1
            out_stream.write(readline_to_write(m_read, paired))
            m_read = update_m()
    elif bool(p_read):
        while bool(p_read):
            pat_only += 1
            out_stream.write(readline_to_write(p_read, paired))
            p_read = update_p()                

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
    main()
