#! /usr/bin/env python

import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import gzip
import argparse
import re
import copy


def parseArgs():
    '''
    Parses arguments. An input CDS file and an output file are required. 
    An output file is optional. If an output file is not specified, output will be writted to standard out.
    '''
    import os
    parser = argparse.ArgumentParser(description='Find Orthologs in Two Files.')
    parser.add_argument("-c",help="Input CDS file",action="store", dest="cdsFile", required=True)
    #parser.add_argument("-v",help="Input VCF file",action="store", dest="vcfFile", required=True)
    #parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
    parser.add_argument("-i",help="Input VCF File",action="store", dest="input", required=False)
    parser.add_argument("-o",help="Output Directory",action="store",dest="output", required=True)
    #parser.add_argument("-t",help="Number of Cores",action="store",dest="threads",type=int, default=-1, required=False)
    args = parser.parse_args()
    if not os.path.isfile(args.cdsFile):
        print (args.cdsFile, "is not a correct file path!")
        sys.exit()
    #if not os.path.isfile(args.vcfFile):
    #    print args.vcfFile, "is not a correct file path!"
    #    sys.exit()

    return args

def readCDS(inputFile):
    cds = dict()
    curHeader = 0
    headerDict =dict() #so multiple header lines can be identical
    headerPos = dict()
    totalSeq = ""
    forward= True
    positions = []
    chrom = ""
    for line in open(inputFile):
        if line[0] == '>':
            if totalSeq != "":
                if not forward:
                    nucSeq = str(Seq(totalSeq).reverse_complement())
                    totalSeq = nucSeq[::-1]
                for x in range(len(positions)):
                    if not positions[x] in cds[chrom]:
                        cds[chrom][positions[x]] = totalSeq[x].upper()
                    if forward:
                        headerPos[chrom][curHeader]=tuple(["+",positions])
                    else:
                        headerPos[chrom][curHeader] =tuple(["-",positions])
            headerDict[curHeader] = line
            curHeader += 1
            totalSeq = ""
            forward = True
            positions = []
            continue
        info = line.strip().split("|")
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        seq = info[3]
        totalSeq +=seq
        if not chrom in headerPos:
            headerPos[chrom] = dict()
        if not (curHeader) in headerPos[chrom]:
            headerPos[chrom][curHeader] = []
        if not chrom in cds:
            cds[chrom] = dict()
        if min(start,end) == start:
            for x in range(start,end+1,1):
                positions.append(x)
        else:
            for x in range(start,end-1,-1):
                forward = False
                positions.append(x)
    if not forward:
        nucSeq = str(Seq(totalSeq).reverse_complement())
        totalSeq = nucSeq[::-1]
    for x in range(len(positions)):
        if not positions[x] in cds[chrom]:
            cds[chrom][positions[x]] = totalSeq[x].upper()
        if forward:
            headerPos[chrom][curHeader]=tuple(["+",positions])
        else:
            headerPos[chrom][curHeader]=tuple(["-",positions])
    return headerPos,cds,headerDict


def readVCF(vcf):
    if vcf.endswith(".gz"):
        vcfFile = gzip.open(vcf,'r')
    else:
        vcfFile = open(vcf,'r')
    allAlts = {}
    for line in vcfFile:
        if isinstance(line,bytes):
            line=line.decode('UTF-8')
        if line[0] == '#':
            if line[1]=='#':
                continue
            ###indNums =line.strip().split('\t')[9:]
            continue

        info = line.strip().split("\t")
        chrom = info[0]
        pos = int(info[1])
        ref = info[3]
        #alt = [ref]
        #alt.extend(info[4].split(",")) #all alt's as a list
        alt=info[4].split(",") #all alt's as a list
        alt = tuple(alt)
        if not chrom in allAlts:
            allAlts[chrom] = dict()
        allAlts[chrom][pos] = tuple([alt, "\t".join(info[0:8]),ref])

    vcfFile.close()
    return allAlts

def writeFile(allAlts,cds,headerPos,headerDict,output):
    outputFile = open(output,'w')
    for chrom, headers in headerPos.items():
        for header,tuplePos in headers.items():
            strand = tuplePos[0]
            numList = tuplePos[1]
            for num in numList:
                if chrom in allAlts and num in allAlts[chrom]:
                    outputFile.write(headerDict[header])
                    if strand == "+":
                         for curPos in numList:
                            outputFile.write(cds[chrom][curPos])
                    else:
                        for curPos in numList:
                            nuc = Seq(cds[chrom][curPos])
                            outputFile.write(str(nuc.reverse_complement()))
                    outputFile.write("\n")
                    for alternate in allAlts[chrom][num][0]:
                        outputFile.write(">" +allAlts[chrom][num][1] +"\t" + alternate + "\n")
                        ref = allAlts[chrom][num][2]
                        removedNums = set()
                        if len(ref) >1:
                            for x in range(1,len(ref)):
                                removedNums.add(num +x)
                        if strand == "+":
                            for curPos in numList:
                                if curPos in removedNums:
                                    pass
                                elif curPos!= num:
                                    outputFile.write(cds[chrom][curPos])
                                else:
                                    if "<CN" in alternate:
                                        copy_num = int(alternate.split("<CN")[1].replace(">",""))
                                        outputFile.write(copy_num*ref)
                                    else:
                                        outputFile.write(alternate)
                        else:
                            for curPos in numList:
                                if curPos in removedNums:
                                    pass
                                elif curPos!= num:
                                    nuc = Seq(cds[chrom][curPos])
                                    outputFile.write(str(nuc.reverse_complement()))
                                else:
                                    if "<CN" in alternate:
                                        copy_num = int(alternate.split("<CN")[1].replace(">",""))
                                        nuc = Seq(copy_num*ref)
                                        outputFile.write(str(nuc.reverse_complement()))
                                    else:
                                        nuc = Seq(alternate)
                                        outputFile.write(str(nuc.reverse_complement()))
                        outputFile.write("\n")


if __name__ =='__main__':
    '''
    Main
    '''
    #freeze_support()
    args = parseArgs()
    headerPos,cds,headerDict = readCDS(args.cdsFile)
    print ("Read CDS")
    allAlts = readVCF(args.input)
    print ("Read VCF")
    writeFile(allAlts,cds,headerPos,headerDict,args.output)



