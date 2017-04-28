'''
"--GCMat" -- path to gene counts matrix
"--gtf" -- path to gtf file
"--trIDs" -- path to file with column of transcripts
"--o" -- path to folder for output, with no '/' at the end
'''

import argparse
import re 


def ObtainGCDict(GCMat):
    res = {}
    with open (GCMat, 'r') as f:
        for line in f:
            data = line.strip().split()
            res[data[0]] = data[1]
    return res


def ObtainTGDict(trIDs, gtf):
    res = {}
    with open (trIDs, 'r') as f:
        with open (gtf, 'r') as g:
            gtf_line = g.read()
            for line in f:
                tr = line.strip().split()[0]
                if tr not in res: 
                    pattern = 'gene_id "[a-zA-Z0-9_]*"; transcript_id "' + tr + '"'
                    gen = re.search(pattern, gtf_line).group(0).split('"')[1]
                    res[tr] = gen
    return res


def CountsIDsGtoTallinone(GCMat, trIDs, gtf, ofolder):
    GenCounts = ObtainGCDict(GCMat)
    TrGen = ObtainTGDict(trIDs, gtf)
    
    with open(GCMat.split('.')[0] + '_TRANSCRIPT' + '.txt', 'w') as W:
        for item in TrGen:
            W.write(item + '\t' + GenCounts[TrGen[item]] + '\n')

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--GCMat", required=True, help="path to gene counts matrix")
    parser.add_argument("--gtf", required=True, help="path to gtf file")
    parser.add_argument("--trIDs", required=True, help="path to file with column of transcripts")
    parser.add_argument("--o", required=True, help="path to folder for output, with no '/' at the end")
    args = parser.parse_args()

    CountsIDsGtoTallinone(args.GCMat, args.trIDs, args.gtf, args.o)
    
if __name__ == "__main__":
    main()
