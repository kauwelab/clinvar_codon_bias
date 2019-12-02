#! /usr/bin/env python
import sys
import argparse
import re
from collections import Counter

def getCodonAversion(sequence, possibleCodons):
    '''
    Takes two arguments: A DNA or RNA sequence, and a set of all possible codons.
    Subtracts a set of all codons in the sequence from the set of all possible codons.
    Returns a tuple of the subtracted set. This tuple represents all codons not found in the sequence.
    '''
    foundCodons= set(re.findall("...",sequence))
    motif = tuple(sorted(list(possibleCodons - foundCodons)))
    return tuple(motif)

def makeAllPossibleCodons(rna,co_trna):
    '''
    Input: rna is a flag to specify if the sequence is DNA or RNA. co_trna is a flag to
        specify co_tRNA codon pairing (i.e., same amino acid formed)
    Returns a set of all 64 possible codons (DNA or RNA) or all 20 amino acids.
    '''
    if co_trna:
        return set(['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
    from itertools import product
    codons = product("ACGT",repeat=3)
    if rna:
        codons = product("ACGU",repeat=3)
    codonsComb = set()
    for c in codons:
        codonsComb.add("".join(c))
    return codonsComb
def parseArgs():
    '''
    Argument parsing is done.
    Required to have an input file.
    '''
    parser = argparse.ArgumentParser(description='Find Differences in codon pairing and codon aversion')
    #parser.add_argument("-t",help="Number of Cores",action="store",dest="threads",default=0,type=int, required=False)
    parser.add_argument("-i",help="Input VCF File",action="store", dest="inputVCF", required=False)
    parser.add_argument("-a",help="Output Codon Aversion File",action="store",dest="output_aversion", required=False)
    parser.add_argument("-o",help="Output Codon Pairing File",action="store",dest="output_pairing", required=False)
    parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
    parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
    parser.add_argument("-b",help="Both Identical and Co-tRNA codon pairing",action="store_true",dest="both", required=False)
    parser.add_argument("-e",help="Either Identical or Co-tRNA codon pairing",action="store_true",dest="either", required=False)
    parser.add_argument("-rna",help="Flag for RNA sequences",action="store_true",dest="rna", required=False)
    parser.add_argument("-l",type=str, help="Codon Table. Default: Standard",action="store",dest="codon_table", default="Standard", required=False)
    args = parser.parse_args()

    if args.co_trna and args.both:
        sys.stdout.write("You cannot use both the co_trna (-c) and both (-b) flags.\n")
        sys.exit()

    return args

def getPairs(seq,orderedCodons):
    footprint = args.footprint
    pairs = dict()
    codons = []
    dna_codons = []
    if args.both or args.co_trna:
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        from Bio.Alphabet import generic_rna
        if args.rna:
            rna = Seq(seq,generic_rna)
            aa = str(rna.translate(table=args.codon_table))
            codons = re.findall(".",aa)
        else:
            sequence = Seq(seq,generic_dna)
            aa = str(sequence.translate(table=args.codon_table))
            codons = re.findall(".",aa)
    else:
        codons = re.findall("...",seq)
    if args.co_trna: #To ensure that identical codon pairing does not form the amino acid
            dna_codons = re.findall("...",seq)

    lastFound = dict() #key= codon, value= position of last found codon with pairing #For co-trna: key = codon (amino acid), value= dict() where key=dna_codon (codon) and value = last position of it
    foundPairing = Counter()
    for x in range(len(codons)):
        curCodon = codons[x]
        if not curCodon in orderedCodons:
            continue
        if not args.co_trna:
            if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint): #Must be >= because if footprint is 2 and AAA is found at positions 3 and 4, 4-3 =1, which is 1 less than footprint size.
                lastFound[curCodon] =x
                continue
        else:
            if not curCodon in lastFound:  
                lastFound[curCodon] = dict()
                lastFound[curCodon][dna_codons[x]] =x
                continue
            closestPos = -100
            for key,value in lastFound[curCodon].items():
                if key == dna_codons[x]:
                    continue
                if value >closestPos:
                    closestPos = value
            if (x - closestPos) >= footprint:
                lastFound[curCodon][dna_codons[x]] =x
                continue
        foundPairing[curCodon] +=1
        if args.co_trna:
            lastFound[curCodon][dna_codons[x]] =x
            continue
        lastFound[curCodon] = x
    return foundPairing

if __name__ =='__main__':
    '''
    Main.
    '''
    args = parseArgs()
    inputVCF = open(args.inputVCF)
    outputCodonPairing = open(args.output_pairing,'w')
    outputCodonAversion = open(args.output_aversion,'w')
    
    lastHead = ""
    lastSeq = ""
    lastCodonAversion = set()
    lastCodonPairing = set()
    head = inputVCF.readline()
    orderedAA = makeAllPossibleCodons(args.rna, (args.co_trna | args.both | args.either))
    orderedCodons = makeAllPossibleCodons(False, False) # for codon aversion
    count =0
    isoformChanges = dict()
    types= ['CLNSIG=Benign','CLNSIG=Likely_benign','CLNSIG=Likely_pathogenic','CLNSIG=Pathogenic','CLNSIG=risk_factor','CLNSIG=not_provided','CLNSIG=Uncertain_significance','CLNSIG=drug_response']
    support=['CLNREVSTAT=no_assertion','CLNREVSTAT=criteria_provided,_single_submitter','CLNREVSTAT=criteria_provided,_multiple_submitters','CLNREVSTAT=reviewed_by_expert_panel']
    for s in support:
        isoformChanges[s] = dict()
        for t in types:
            isoformChanges[s][t] =dict()

    while head !="":
        notWrite = False
        seq = inputVCF.readline().strip()
        if head[1:4] == 'ID=':
            lastHead = head
            lastSeq = seq
            if len(seq) %3!=0 or "." in seq:
                notWrite=True
                head = inputVCF.readline()
                count +=1
                continue
            if args.either:
                args.co_trna=True
                lastCodonPairing = getPairs(seq,orderedAA)
                args.co_trna=False
                lastCodonPairing |= getPairs(seq,orderedCodons)
            else:
                lastCodonPairing = getPairs(seq,orderedAA)
            #lastCodonPairing = getPairs(seq,orderedAA)
            lastCodonAversion = getCodonAversion(seq, orderedCodons)
        else:
            for typeVar in types:
                if typeVar in head:
                    t = typeVar
                    break
            for supportVar in support:
                if supportVar in head:
                    s = supportVar
                    break
            if notWrite:
                notWrite=False
                head = inputVCF.readline()
                continue
            if len(seq) %3 !=0 or "." in seq:
                count+=1
                head = inputVCF.readline()
                continue
            if t =='':
                print(head)
            if s =='':
                print(head)
            if not head in isoformChanges[s][t]:
                isoformChanges[s][t][head] =dict()
                isoformChanges[s][t][head]['aversion'] =[0,0,0] #less,same,more
                isoformChanges[s][t][head]['pairing'] =[0,0,0] #less,same,more
                isoformChanges[s][t][head]['both'] =0 #count
            if args.either:
                args.co_trna=True
                codonPairing = getPairs(seq,orderedAA)
                args.co_trna=False
                codonPairing |= getPairs(seq,orderedCodons)
            else:
                codonPairing = getPairs(seq,orderedAA)
            codonAversion = getCodonAversion(seq, orderedCodons)
            changedPairing =False
            if lastCodonPairing != codonPairing:
                changedPairing =True
                difference =sum(codonPairing.values())-sum(lastCodonPairing.values())
                if difference >1: # possible for codon to pair multiple times
                    difference =1
                if difference <-1:
                    difference = -1
                isoformChanges[s][t][head]['pairing'][difference+1] +=1
                outputCodonPairing.write(lastHead +str(lastCodonPairing) + "\n" + head +str(codonPairing) +"\n")
            if lastCodonAversion != codonAversion:
                if changedPairing:
                    isoformChanges[s][t][head]['both']+=1

                difference = len(lastCodonAversion)-len(codonAversion)
                isoformChanges[s][t][head]['aversion'][difference+1] +=1
                outputCodonAversion.write(lastHead +str(lastCodonAversion) + "\n" + head +str(codonAversion) +"\n")
        head = inputVCF.readline()
    for s in isoformChanges:
        for t in isoformChanges[s]:
            isoformChanges[s][t]["aversion"] =[0,0,0,0,0] #less,same,more,inconclusive,all
            isoformChanges[s][t]["pairing"] =[0,0,0,0,0]
            isoformChanges[s][t]["totalAversion"] =[0,0,0] #less,same,more,inconclusive
            isoformChanges[s][t]["totalPairing"] =[0,0,0]
            isoformChanges[s][t]["both"] =[0,0] #totalIsoformChanges, totalCombined
            for head in isoformChanges[s][t]:
                if head in ['aversion','pairing','totalAversion','totalPairing','both']:
                    continue
                numBoth = isoformChanges[s][t][head]['both']
                if numBoth >0:
                    isoformChanges[s][t]['both'][0] +=numBoth
                    isoformChanges[s][t]['both'][1] +=1
                for x in ['aversion','pairing']:
                    num = Counter(isoformChanges[s][t][head][x]).most_common()
                    del num[0]
                    isoformChanges[s][t][x][4] +=1 #all synonymous variants
                    if len(num) == 0:
                        continue
                    name = "total" + x.capitalize()
                    isoformChanges[s][t][name][0] +=isoformChanges[s][t][head][x][0]
                    isoformChanges[s][t][name][1] +=isoformChanges[s][t][head][x][1]
                    isoformChanges[s][t][name][2] +=isoformChanges[s][t][head][x][2]
                    num = num[0]
                    if num[1] !=1:
                        isoformChanges[s][t][x][3] +=1
                        continue
                    index = isoformChanges[s][t][head][x].index(num[0])
                    isoformChanges[s][t][x][index] +=1

    for x in ['aversion','pairing']:
        print(x)
        for s in isoformChanges:
            print(s)
            for t in isoformChanges[s]:
                name = "total" + x.capitalize()
                isoformChanges[s][t][name][0]
                print(",".join(list(map(str,[t,isoformChanges[s][t][x][4],sum(isoformChanges[s][t][x][0:4]),isoformChanges[s][t][x][0],isoformChanges[s][t][x][1],isoformChanges[s][t][x][2],isoformChanges[s][t][x][3],sum(isoformChanges[s][t][name]),isoformChanges[s][t][name][0],isoformChanges[s][t][name][1],isoformChanges[s][t][name][2]])))) #totalSynonymous,numX,less,same,more,inconclusive,totalIsoform,isoformLess,isoformSame,isofomMore
            print ("")
        print ("")
    for x in ['both']:
        print(x)
        for s in isoformChanges:
            print(s)
            for t in isoformChanges[s]:
                name = "both"
                isoformChanges[s][t][name][0]
                print(",".join(list(map(str,[t,isoformChanges[s][t][x][0],isoformChanges[s][t][x][1]]))))#Both num isoforms, Both num total
            print ("")
        print ("")

                



