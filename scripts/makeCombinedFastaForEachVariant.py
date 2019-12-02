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
    curHeader = 1
    headerDict =dict()
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
                        headerPos[chrom][curHeader-1]=tuple(["+",positions])
                    else:
                        headerPos[chrom][curHeader-1] =tuple(["-",positions])
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
        if not (curHeader-1) in headerPos[chrom]:
            headerPos[chrom][curHeader-1] = []
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
            headerPos[chrom][curHeader-1]=tuple(["+",positions])
        else:
            headerPos[chrom][curHeader-1]=tuple(["-",positions])
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
        #ref = info[3]
        #alt = [ref]
        #alt.extend(info[4].split(",")) #all alt's as a list
        alt=info[4].split(",") #all alt's as a list
        alt = tuple(alt)
        if not chrom in allAlts:
            allAlts[chrom] = dict()
        allAlts[chrom][pos] = tuple([alt, "\t".join(info[0:8])])

    vcfFile.close()
    return allAlts

def writeFile(allAlts,cds,headerPos,headerDict,output):
    outputFile = open(output,'w')
    for chrom, headers in headerPos.items():
        for header,tuplePos in headers.items():
            refSeq = ""
            refPos = {}
            strand = tuplePos[0]
            numList = tuplePos[1]
            for num in numList:
                if chrom in allAlts and num in allAlts[chrom]:
                    completelySynonymous = True
                    toWrite = ""
                    for alternate in allAlts[chrom][num][0]:
                        if refSeq=="":
                            count = 0
                            #outputFile.write(headerDict[header])
                            if strand == "+":
                                 for curPos in numList:
                                    refPos[curPos] = count
                                    count +=1
                                    #outputFile.write(cds[chrom][curPos])
                                    refSeq += cds[chrom][curPos]
                            else:
                                originalSeq = ""
                                for curPos in numList:
                                    refPos[curPos] = count
                                    count +=1
                                    originalSeq +=cds[chrom][curPos]
                                    #nuc = Seq(cds[chrom][curPos])
                                    #outputFile.write(str(nuc.reverse_complement()))
                                nuc = Seq(originalSeq)
                                refSeq = str(nuc.reverse_complement())[::-1]
                        #outputFile.write("\n")
                        #outputFile.write(">" +allAlts[chrom][num][1] +"\t" + alternate + "\n")
                        altSeq = ""
                        if strand == "+":
                            pos = refPos[num]
                            altSeq = refSeq[0:pos] + alternate + refSeq[pos+1:]
                            #for curPos in numList:
                            #    elif curPos!= num:
                            #        #outputFile.write(cds[chrom][curPos])
                            #        altSeq += cds[chrom][curPos]
                            #    else:
                            #        altSeq += alternate
                        else:
                            pos = refPos[num]
                            altSeq = refSeq[0:pos] + str(Seq(alternate).reverse_complement()) + refSeq[pos+1:]
                            #altSeq[num] = str(Seq(alternate).reverse_complement())
                            #for curPos in numList:
                            #    elif curPos!= num:
                            #        nuc = Seq(cds[chrom][curPos])
                            #        #outputFile.write(str(nuc.reverse_complement()))
                            #        altSeq += str(nuc.reverse_complement())
                            #    else:
                            #        nuc = Seq(alternate)
                            #        altSeq +=str(nuc.reverse_complement())
                            #        #outputFile.write(str(nuc.reverse_complement()))
                        refProt = str(Seq(refSeq).translate())
                        altProt = str(Seq(altSeq).translate())
                        if refProt!=altProt:
                            sys.stdout.write(">" +allAlts[chrom][num][1] +"\t" + alternate + "\n")
                            completelySynonymous = False
                            break
                        toWrite += headerDict[header]
                        toWrite += (refSeq +"\n")
                        toWrite += (">" +allAlts[chrom][num][1] +"\t" + alternate + "\n")
                        toWrite += (altSeq +"\n")
                    if completelySynonymous:
                        outputFile.write(toWrite)


                        ###if refProt == altProt: #make sure it is synonymous for the isoform
                        ###    outputFile.write(headerDict[header])
                        ###    outputFile.write(refSeq +"\n")
                        ###    outputFile.write(">" +allAlts[chrom][num][1] +"\t" + alternate + "\n")
                        ###    outputFile.write(altSeq +"\n")
                            
                        #outputFile.write("\n")


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



