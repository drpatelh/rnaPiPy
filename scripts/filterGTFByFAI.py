#!/usr/bin/env python

import argparse
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import utilityFunctions as uF

############################################
## PARSE ARGUMENTS                        ##
############################################

Description = 'Filter GTF file to include chromosomes in FAI file.'
Epilog = """Example usage: python filterGTFByFAI.py <FAI_IN> <GTF_IN> <GTF_OUT>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('FAI_IN', help="Full path to FAI file.")
argParser.add_argument('GTF_IN', help="Full path to GTF input file.")
argParser.add_argument('GTF_OUT', help="Full path to GTF output file.")
args = argParser.parse_args()

############################################
## FILTER GTF FILE                        ##
############################################

if exists(args.FAI_IN):
    
    ChromIDs = []    
    fin = open(args.FAI_IN,'r')
    for line in fin.readlines():
        lspl = [x.strip() for x in line.strip().split('\t')]
        ChromIDs.append(lspl[0])
    fin.close()
    
    if exists(args.GTF_IN):
                
        uF.makedir(dirname(args.GTF_OUT))
        fin = open(args.GTF_IN,'r')
        fout = open(args.GTF_OUT,'w')
        while True:
            line = fin.readline()
            if line:
                if not line[0] == '#':
                    lspl = line.strip().split('\t')
                    if lspl[0] in ChromIDs:
                        fout.write(line)
                else:
                    fout.write(line)
            else:
                fin.close()
                fout.close()
                break

############################################
############################################
############################################
############################################