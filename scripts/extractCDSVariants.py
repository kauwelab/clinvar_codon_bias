#! /usr/bin/env python3
import sys
import gzip
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description='Find Orthologs in Two Files.')
    parser.add_argument("-i",help="Input CDS file",action="store", dest="cds", required=True)
    parser.add_argument("-c",help="Input Gzipped ClinVar File",action="store",dest="clinvar", required=True)
    parser.add_argument("-o",help="Output File",action="store",dest="output", required=True)
    args = parser.parse_args()
    return args


def readCDSFile(CDS_file):
    input1 = open(CDS_file,'r') #CDS file
    everyNum = dict()
    allSeq = dict()
    for line in input1:
        if line[0] == '>':
            continue
        info = line.strip().split("|")
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        seq = info[3]
        if not chrom in everyNum:
            #everyNum[chrom] = dict()
            everyNum[chrom] = set()
        if min(start,end) == start:
            for x in range(start,end+1,1):
                everyNum[chrom].add(x)
                #everyNum[chrom][x] = seq[x-start]
        else:
            for x in range(start,end-1,-1):
                everyNum[chrom].add(x)
                #everyNum[chrom][x] = seq[len(seq)-(x-start)-1]
    input1.close()
    return everyNum

def findSynonymousVariants(clinvar_file,output_file,everyNum):
    if clinvar_file.endswith(".gz") or clinvar_file.endswith(".gzip"):
        clinvarFile = gzip.open(clinvar_file,'r')
    else:
        clinvarFile = open(clinvar_file,'r')
    output = open(output_file,'w')
    for line in clinvarFile:
        if not isinstance(line,str):
            line = line.decode()
        if line[0] =='#':
            continue
        #if "intron_variant" in line:
        #    continue
        #if onlySynonymous:
        #    if not "|synonymous_variant;" in line:
        #        continue
        info = line.strip().split("\t")
        chrom = info[0] 
        pos = int(info[1])
        ref = info[3]
        alt = info[4]
        if len(alt) >1 or len(ref) >1:
            continue
        if alt == ".":
            continue
        if chrom in everyNum:
            if pos in everyNum[chrom]:
                output.write(line)

    clinvarFile.close()
    output.close()
    
if __name__ =='__main__':
    '''
    Main
    '''
    args = parseArgs()
    everyNum = readCDSFile(args.cds)
    findSynonymousVariants(args.clinvar,args.output,everyNum)

