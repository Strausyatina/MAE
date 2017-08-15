# input: fasta_to_mutate, mutation_rate(one_per_what?), fasta_output
# output : fasta with randomly inserted mutations,
#	   list of triplets (pos, old, new) 
# usage: python3 FastaRandomMutations.py --fasta smth.fa --ofasta smth_mut.fa --mutrate 1000
# note: ORIENTED FOR USAGE ON (big seqs) CHROMOSOMES


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import numpy as np
import argparse

replacements = {'A':['T','C','G'], 'T':['A','C','G'], 'C':['A','T','G'], 'G':['A','T','C'], 'N':['N']}

def mutateRecord (record, mutRate):
    l = len(record)
    n = int(math.ceil(l/int(mutRate)))
    seq = str(record.seq)
    subseq = np.random.choice(l-1, n, replace=False)
    subseq.sort()
    
    print(l, n, subseq)

    mutsubseq = [(i, seq[i], np.random.choice(replacements[seq[i]], 1)[0]) for i in subseq]
    newseq = seq[:mutsubseq[0][0]]
    
    #print(mutsubseq)

    for i in range(n):
         if (i != n-1):
            newseq = newseq + mutsubseq[i][2] + seq[(mutsubseq[i][0]+1):mutsubseq[i+1][0]]
         else:
            newseq = newseq + mutsubseq[i][2] + seq[(mutsubseq[i][0]+1):]
         if (i % 1000 == 0):
            print("another 1000")

    mutrecord = SeqRecord(Seq(newseq), id = record.id, 
                          name = record.name, description = record.description)
    return (mutrecord, mutsubseq)

def mutateFasta (inputFa, outputFa, mutRate):
    with open (outputFa, 'w') as outfa:
        with open (outputFa+".log", 'w') as outlog:
            for record in SeqIO.parse(inputFa, "fasta"):
                mutrecord , mutsubseq = mutateRecord(record, mutRate)
                SeqIO.write(mutrecord, outfa, "fasta")
                
                outlog.write(mutrecord.description + '\n')
                for triplet in mutsubseq:
                    outlog.write('\t'.join([str(item) for item in triplet]) + '\n')

def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="path to fasta for mutating")
    parser.add_argument("--ofasta", required=True, help="path to mut_fasta output")
    parser.add_argument("--mutrate", required=True, help="rate: one mutation per _what_")
    args = parser.parse_args()
    
    mutateFasta(args.fasta, args.ofasta, args.mutrate)
    
if __name__ == "__main__":
    main()

