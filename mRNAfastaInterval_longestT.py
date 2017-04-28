'''
"--fasta" -- path to fasta
"--gtf" -- path to gtf file
"--chrom" -- chromosome
"--s","--e" -- interval on chromosome
"--ofasta" -- path to subfasta outbut
"--oGT" -- path to file with pairs Gene '\t' Longest Transcript
'''

import argparse

def getGenesOnSegment (gtf, chrom, I):
    geneSet = set()
    with open (gtf, 'r') as S:
        for line in S:
            if line.startswith('#'):
                continue
            row = line.strip().split('\t')
            if int(row[4])<=I[0] or I[1]<=int(row[3]):
                    continue
            if row[0] == str(chrom) and row[2] == "exon":
                infoList = [item.split('"') for item in row[8].replace(' ','').split(';')[:-1]]
                info = {item[0]:item[1] for item in infoList}
                if 'gene_id' in info.keys():
                    geneSet.add(info['gene_id'])
    return list(geneSet)

def obtainGenesAndLongestTranscripts (gtf, genes, GToutput): 
    with open (gtf, 'r') as S:
        coords = [10**20,0]
        dictTrLen = {}
        dictTrCoords = {}
        dictG = {}
        for line in S:
            if line.startswith('#'):
                continue
            row = line.strip().split('\t')
            if row[2] == "exon":
                infoList = [item.split('"') for item in row[8].replace(' ','').split(';')[:-1]]
                info = {item[0]:item[1] for item in infoList}
                if 'gene_id' in info.keys():
                    geneID = info['gene_id']
                    if geneID in genes:
                        s , e = int(row[3]) , int(row[4])
                        if s < coords[0]:
                            coords[0] = s
                        if e > coords[1]:
                            coords[1] = e
                        #print(geneName, ",", info["gene_id"] ,",", 
                        #      info["transcript_id"], " : ", 
                        #      s, "-", e)
                        if info["gene_id"] in dictG.keys():
                            if info["transcript_id"] not in dictG[info["gene_id"]]:
                                dictG[info["gene_id"]].append(info["transcript_id"])
                        else: 
                            dictG[info["gene_id"]] = [info["transcript_id"]]
                        if info["transcript_id"] in dictTrLen.keys():
                            dictTrLen[info["transcript_id"]] += e-s+1
                            dictTrCoords[info["transcript_id"]].append((row[0],s,e))
                        else:
                            dictTrLen[info["transcript_id"]] = e-s+1
                            dictTrCoords[info["transcript_id"]] = [(row[0],s,e)]
    #print(coords)        
    #129: [4837602, 5812851]
    #CAST: [4782800, 5745899]
    #38: [7811054, 8708742]
    
    GT = {}
    GTcoords = {}
    with open (GToutput, 'w') as f:
        f.write ("#GeneID" + '\t' + "LongestTranscriptID" + '\n')
        for item in dictG:
            maxL = 0
            maxT = ''
            for jtem in dictG[item]:
                l = int(dictTrLen[jtem])
                if maxL < l:
                    maxL , maxT , maxTcoords = l , jtem , dictTrCoords[jtem]
                #print(item, jtem, dictTrLen[jtem])
            GT[item] = maxT
            GTcoords[maxT] = maxTcoords
            f.write(item + '\t' + maxT + '\n')
    #print(GT)        
    return GT , GTcoords

def getGTFsubGeneEXON (gtf, output_name, GToutput, genes): 
    GenTransDict = obtainGenesAndLongestTranscripts (gtf, genes, GToutput)[0]
    #print(obtainGenesAndLongestTranscripts (gtf, genes, GToutput)[1])
    #print(GenTransDict)
    with open (gtf, 'r') as S:
        with open (output_name, 'w') as O:
            for line in S:
                if line.startswith('#'):
                    continue
                row = line.strip().split('\t')
                #if row[0] == '15' and row[2] == "exon":
                if row[2] == "exon":
                    infoList = [item.split('"') for item in row[8].replace(' ','').split(';')[:-1]]
                    info = {item[0]:item[1] for item in infoList}
                    if 'gene_id' in info.keys():
                        gID = info['gene_id']
                        if gID in GenTransDict.keys() and info["transcript_id"] == GenTransDict[gID]:
                            #print(line)
                            O.write(line)


def getGTFsubGeneEXONonSegm (gtf, output_name, GToutput, chrom, I):
    getGTFsubGeneEXON (gtf, output_name, GToutput, getGenesOnSegment(gtf,chrom,I))

def StringSplitN (string, N):
    newstr = ''
    counter = 0  
    for c in string:
        if counter % N == 0:
            newstr += '\n'
        newstr += c
        counter += 1
    return newstr

def GetTranscriptFasta (gID, trID, coords, fasta):
    '''
    Takes coordinates and returns subFASTA
    '''
    if not coords:
        print('no intervals provided')
        return 
    
    chrom = coords[0][0]
    coords = sorted(coords)
    
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
            
    name = '>' + gID + ':' + trID
    return name , transcript

def GetIntervalMRNAfasta(fasta, gtf, chrom, I, GLToutFile, subFastaFile):
    GeneLongTrDict , TrCoords = obtainGenesAndLongestTranscripts(gtf,getGenesOnSegment(gtf,chrom,I),GLToutFile)
    
    with open (subFastaFile, 'w') as OUT:
        for gID in GeneLongTrDict:
            trID = GeneLongTrDict[gID]
            coords = TrCoords[trID]
            name , tr = GetTranscriptFasta (gID, trID, coords, fasta)
            OUT.write(name + StringSplitN(tr, 60) + '\n')            
    print("fasta subfile is written")

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="path to fasta")
    parser.add_argument("--gtf", required=True, help="path to gtf file")
    parser.add_argument("--chrom", required=True, help="cromosome")
    parser.add_argument("--s", required=True, help="start position")
    parser.add_argument("--e", required=True, help="end position")
    parser.add_argument("--ofasta", required=True, help="path to subfasta outbut")
    parser.add_argument("--oGT", required=True, help="path to file with pairs Gene '\t' Longest Transcript")
    args = parser.parse_args()
    
    GetIntervalMRNAfasta(args.fasta, args.gtf, args.chrom, [int(args.s),int(args.e)], args.oGT, args.ofasta)
    
if __name__ == "__main__":
    main()

