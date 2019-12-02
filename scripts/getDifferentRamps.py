#! /usr/bin/env python

import sys
import argparse
def parseArgs():
    parser = argparse.ArgumentParser(description='Find Variants that destroy or create a ramp sequence.')
    parser.add_argument("-i",help="Input Ramps file",action="store", dest="input", required=True)
    parser.add_argument("-r",help="Input Original File",action="store",dest="original", required=True)
    parser.add_argument("-o",help="Output File",action="store",dest="output", required=False)
    args = parser.parse_args()
    return args

def readRampFile(args):
    rampFile = open(args.input) 
    orderedRamps = []
    head = rampFile.readline()
    while head != "":
        seq = rampFile.readline()
        orderedRamps.append(head)
        head = rampFile.readline()
    rampFile.close()
    return orderedRamps

def readOriginalFile(args):
    originalOrder = open(args.original,'rb')
    originalHead = originalOrder.readline().decode()
    varOriginal = {}
    originalVar = {}
    while originalHead != "":
        originalOrder.readline()
        varHead = originalOrder.readline().decode()
        originalOrder.readline()
        if not varHead in varOriginal:
            varOriginal[varHead] = []
        varOriginal[varHead].append(originalHead)
        if not originalHead in originalVar:
            originalVar[originalHead] = []
        originalVar[originalHead].append(varHead)
        originalHead = originalOrder.readline().decode()
    originalOrder.close()
    return varOriginal,originalVar

def getChangeRampExistence(orderedRamps,varOriginal,originalVar):
    createOrDestroy =set()
    for ramp in orderedRamps:
        if ramp.startswith(">ID="):
            ref = ramp
            posVars = originalVar[ref]
            for var in posVars:
                if not var in orderedRamps:
                    createOrDestroy.add(var)
                    #output.write(var)
                    continue
        else:
            var = ramp
            refs = varOriginal[var]
            for ref in refs:
                if not ref in orderedRamps:
                    createOrDestroy.add(var)
                    #output.write(var)
    return createOrDestroy

    

if __name__ =='__main__':
    '''
    Main
    '''
    args= parseArgs()
    orderedRamps = readRampFile(args)
    varOriginal,originalVar = readOriginalFile(args)
    createOrDestroy = getChangeRampExistence(orderedRamps,varOriginal,originalVar)
    output = sys.stdout
    if args.output:
        output = open(args.output,'w')
    for var in createOrDestroy:
        output.write(var)

