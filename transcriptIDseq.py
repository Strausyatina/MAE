'''
"--trID" -- desired transcript id
"--fasta" -- path to fasta file
"--gtf" -- path to gtf file
"--o" -- path to folder for output, with no '/' at the end

USAGE:
python transcriptIDseq.py --trID ENSMUST00000082424 --fasta ./GRCm38_68.fa --gtf ./
Mus_musculus.GRCm38.68.gtf --o ./test

Will add to --o folder files: 
1.gtf for transcript
2.fasta of transcript (concatenated exons)
'''


import argparse

def GetTranscriptCoordsGTF (trID, gtf, outFolder):
    '''
    Takes gtf-file and transcript ID;
    Returns gtf for desired transcript and vector of coordinate triplets  
    '''
    oFile = outFolder + '/' + gtf.split('.gtf')[0] + '.' + trID + ".gtf"
    coords = []
    
    with open (gtf, 'r') as GTF:
        with open (oFile, 'w') as OUT:
            for line in GTF:
                if line.startswith('#'):
                    continue
                row = line.strip().split('\t')
                
                if row[2] == "exon":
                    infoList = [item.split('"') for item in row[8].replace(' ','').split(';')[:-1]]
                    info = {item[0]:item[1] for item in infoList}
                    if 'transcript_id' in info.keys():
                        if info['transcript_id'] == trID:
                            coords.append( (row[0], int(row[3]), int(row[4])) )
                            OUT.write(line)
    print("gtf subfile is written")
    return coords

def StringSplitN (string, N):
    newstr = ''
    counter = 0  
    for c in string:
        if counter % N == 0:
            newstr += '\n'
        newstr += c
        counter += 1
    return newstr


def GetTranscriptFasta (trID, coords, fasta, outFolder):
    '''
    Takes coordinates and returns subFASTA
    '''
    oFile = outFolder + '/' + fasta.split('.fa')[0] + '.' + trID + '.fa'
    if not coords:
        print('no intervals provided')
        return 
    
    chrom = coords[0][0]
    coords = sorted(coords)
    print(coords)
    
    transcript = ''
    nR , nI = 0 , 0
    
    with open (fasta, 'r') as FA:
        curchr , curtr = False , False
        for line in FA:
            if nI == len(coords):
                break
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if line.startswith('>'+chrom):
                    print("scanning chromosome")
                    curchr = True
                    continue
                else:
                    curchr = False
                    continue
            if not curchr:
                continue

            # well, you're on right chromosome:
            if curtr:
                if 60*nR < coords[nI][2] <= 60*(nR+1):
                    transcript += line[: (coords[nI][2] - 1) % 60 + 1]
                    curtr = False
                    nI += 1
                else: 
                    transcript += line.strip()
            if not curtr:
                while nI < len(coords) and 60*nR < coords[nI][1] <= 60*(nR+1):
                    if coords[nI][1] < coords[nI][2] <= 60*(nR+1):
                        transcript += line[(coords[nI][1] - 1) % 60 : ((coords[nI][2] - 1) % 60) + 1]
                        nI += 1
                    else:
                        transcript += line.strip()[(coords[nI][1] - 1) % 60:]
                        curtr = True
                        break
            nR += 1
            
    with open (oFile, 'w') as OUT:
        OUT.write('>' + trID + StringSplitN(transcript, 60))            
    print("fasta subfile is written")
    return 

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--trID", required=True, help="desired transcript id")
    parser.add_argument("--fasta", required=True, help="path to fasta file")
    parser.add_argument("--gtf", required=True, help="path to gtf file")
    parser.add_argument("--o", required=True, help="path to folder for output, with no '/' at the end")
    args = parser.parse_args()

    print('gtf format')
    GetTranscriptFasta(args.trID, GetTranscriptCoordsGTF(args.trID, args.gtf, args.o), args.fasta, args.o)
    
    
if __name__ == "__main__":
    main()
