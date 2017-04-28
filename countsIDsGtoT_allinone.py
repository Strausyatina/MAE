'''
"--GCMat" -- path to gene counts matrix
"--GenTrIDs" -- path to file with column of transcripts
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

def CountsIDsGtoTallinone(GCMat, trIDs):
    GenCounts = ObtainGCDict(GCMat)
   
    with open(GCMat.split('.')[0] + '_TRANSCRIPT' + '.txt', 'w') as W:
	with open (trIDs, 'r') as f:
		for line in f:
			if line[0] =='#':
				continue
			cur = line.strip().split()
            		W.write(cur[1] + '\t' + GenCounts[cur[0]] + '\n')

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--GCMat", required=True, help="path to gene counts matrix")
    parser.add_argument("--GenTrIDs", required=True, help="path to file with column of transcripts")
    args = parser.parse_args()

    CountsIDsGtoTallinone(args.GCMat, args.GenTrIDs)
    
if __name__ == "__main__":
    main()
