#! /usr/bin/env python3.7
import sys
import argparse
import os

def parseArgs():
    parser = argparse.ArgumentParser(description='Find Orthologs in Two Files.')
    parser.add_argument("-r",help="Input ramps file",action="store", dest="ramps", required=False)
    parser.add_argument("-a",help="Input codon aversion file",action="store", dest="aversion", required=False)
    parser.add_argument("-p",help="Input codon pairing file",action="store", dest="pairing", required=False)
    parser.add_argument("-o",help="Output Directory",action="store",dest="output", required=True)
    parser.add_argument("-f",help="Only Reference",action="store",dest="reference", required=False)
    args = parser.parse_args()
    if not args.reference:
        if not args.ramps:
            sys.stdout.write("A ramps file is required to combine codon usage biases\n")
            sys.exit()
        else:
            if not os.path.isfile(args.ramps):
                sys.stdout.write("Invalid file path: " + args.ramps + "\n")
                sys.exit()
        if not args.aversion:
            sys.stdout.write("A codon aversion file is required to combine codon usage biases\n")
            sys.exit()
        else:
            if not os.path.isfile(args.aversion):
                sys.stdout.write("Invalid file path: " + args.aversion + "\n")
                sys.exit()
        if not args.pairing:
            sys.stdout.write("A codon pairing file is required to combine codon usage biases\n")
            sys.exit()
        else:
            if not os.path.isfile(args.pairing):
                sys.stdout.write("Invalid file path: " + args.pairing + "\n")
                sys.exit()
    else:
        if not os.path.isfile(args.reference):
            sys.stdout.write("Invalid file path: " + args.reference + "\n")
            sys.exit()


    if not args.output.endswith("/"):
        args.output +="/"
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    return args

def readFile(file_name,reference):
    inputFile = open(file_name)
    types= ['CLNSIG=Benign','CLNSIG=Likely_benign','CLNSIG=Likely_pathogenic','CLNSIG=Pathogenic','CLNSIG=risk_factor','CLNSIG=not_provided','CLNSIG=Uncertain_significance','CLNSIG=drug_response']
    support=['CLNREVSTAT=no_assertion','CLNREVSTAT=criteria_provided,_single_submitter','CLNREVSTAT=criteria_provided,_multiple_submitters','CLNREVSTAT=reviewed_by_expert_panel']
    allVar = dict() #CLNREVSTAT -> CLNSIG -> set(headers)
    for s in support:
        allVar[s] = dict()
        for t in types:
            allVar[s][t] = set()
    for line in inputFile:
        if line.startswith(">ID"):
            continue
        if not line.startswith(">"):
            continue
        for s in support:
            if s in line:
                for t in types:
                    if t in line:
                        allVar[s][t].add(line)
                        break
                break
    inputFile.close()
    return allVar

def getOverlap(allAversion,allPairing,allRamps,outputDir):
    types= ['CLNSIG=Benign','CLNSIG=Likely_benign','CLNSIG=Likely_pathogenic','CLNSIG=Pathogenic','CLNSIG=risk_factor','CLNSIG=not_provided','CLNSIG=Uncertain_significance','CLNSIG=drug_response']
    support=['CLNREVSTAT=no_assertion','CLNREVSTAT=criteria_provided,_single_submitter','CLNREVSTAT=criteria_provided,_multiple_submitters','CLNREVSTAT=reviewed_by_expert_panel']
    for s in support:
        dirName = s.split("=")[1] + "/"
        if "," in dirName:
            dirName = dirName.split(",")[1][1:]
        if not os.path.isdir(outputDir + dirName):
            os.mkdir(outputDir + dirName)
        for t in types:
            fileName = t.split("=")[1] + "/"
            if not os.path.isdir(outputDir + dirName + fileName):
                os.mkdir(outputDir + dirName + fileName)
            aversion_pairing_ramps = allAversion[s][t] & allPairing[s][t] & allRamps[s][t]
            output = open(outputDir + dirName + fileName + "aversion_pairing_ramps",'w')
            for var in aversion_pairing_ramps:
                output.write(var)
            output.close()

            aversion_pairing = (allAversion[s][t] & allPairing[s][t]) - aversion_pairing_ramps
            output = open(outputDir + dirName + fileName + "aversion_pairing",'w')
            for var in aversion_pairing:
                output.write(var)
            output.close()

            aversion_ramps = (allAversion[s][t] & allRamps[s][t]) - aversion_pairing_ramps
            output = open(outputDir + dirName + fileName + "aversion_ramps",'w')
            for var in aversion_ramps:
                output.write(var)
            output.close()

            pairing_ramps = (allPairing[s][t] & allRamps[s][t]) - aversion_pairing_ramps
            output = open(outputDir + dirName + fileName + "pairing_ramps",'w')
            for var in pairing_ramps:
                output.write(var)
            output.close()

            only_pairing = ((allPairing[s][t] - aversion_pairing_ramps) - aversion_pairing) - pairing_ramps
            output = open(outputDir + dirName + fileName + "only_pairing",'w')
            for var in only_pairing:
                output.write(var)
            output.close()

            only_aversion = ((allAversion[s][t] - aversion_pairing_ramps) - aversion_pairing) - aversion_ramps
            output = open(outputDir + dirName + fileName + "only_aversion",'w')
            for var in only_aversion:
                output.write(var)
            output.close()

            only_ramps = ((allRamps[s][t] - aversion_pairing_ramps) - aversion_ramps) - pairing_ramps
            output = open(outputDir + dirName + fileName + "only_ramps",'w')
            for var in only_ramps:
                output.write(var)
            output.close()

def writeReference(ref,outputDir):
    types= ['CLNSIG=Benign','CLNSIG=Likely_benign','CLNSIG=Likely_pathogenic','CLNSIG=Pathogenic','CLNSIG=risk_factor','CLNSIG=not_provided','CLNSIG=Uncertain_significance','CLNSIG=drug_response']
    support=['CLNREVSTAT=no_assertion','CLNREVSTAT=criteria_provided,_single_submitter','CLNREVSTAT=criteria_provided,_multiple_submitters','CLNREVSTAT=reviewed_by_expert_panel']
    for s in support:
        dirName = s.split("=")[1] + "/"
        if "," in dirName:
            dirName = dirName.split(",")[1][1:]
        if not os.path.isdir(outputDir + dirName):
            os.mkdir(outputDir + dirName)
        for t in types:
            fileName = t.split("=")[1] + "/"
            if not os.path.isdir(outputDir + dirName + fileName):
                os.mkdir(outputDir + dirName + fileName)
            output = open(outputDir + dirName + fileName + "reference",'w')
            for var in ref[s][t]:
                output.write(var)
            output.close()


if __name__ =='__main__':
    '''
    Main
    '''
    args = parseArgs()
    if args.reference:
        ref = readFile(args.reference,True)
        writeReference(ref,args.output)
        sys.exit()
    allAversion = readFile(args.aversion,False)
    allPairing=readFile(args.pairing,False)
    allRamps = readFile(args.ramps,False)


    getOverlap(allAversion,allPairing,allRamps,args.output)




