#! /usr/bin/env python
import sys
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='Find Orthologs in Two Files.')
    parser.add_argument("-i",help="Input CDS File",action="store", dest="input", required=True)
    parser.add_argument("-o",help="Output Fasta File",action="store",dest="output", required=True)
    args = parser.parse_args()
    return args

if __name__ =='__main__':
    args= parseArgs()
    inputFasta = open(args.input)
    header = inputFasta.readline()
    everything = {}
    while header != "":
        seq = ""
        gene = header.split("gene=")[1].split(";")[0]
        if not gene in everything:
            everything[gene] = tuple([header,seq])

        nextSeq = inputFasta.readline()
        while nextSeq != "" and nextSeq[0] != ">":
            seq +=nextSeq.strip().split("|")[-1]
            nextSeq = inputFasta.readline()
        seq +="\n"

        if len(seq) > len(everything[gene][1]):
            everything[gene] = tuple([header,seq])
        header = nextSeq
    
    output = open(args.output,'w')
    for header,seq in everything.values():
        output.write(header + seq)
    inputFasta.close()
    output.close()



