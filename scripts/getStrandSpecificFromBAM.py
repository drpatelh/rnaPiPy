#!/usr/bin/env python

import argparse
import ConfigParser
from time import gmtime, strftime
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import utilityFunctions as uF

############################################
## PARSE ARGUMENTS                        ##
############################################

Description = 'Run pysam on BAM file mapped to RSEM TRANSCRIPTOME to work out strand-specificity.'
Epilog = """Example usage: python getStrandSpecificFromBAM.py <CONFIG_FILE> <BAM_FILE> <OUTPREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('BAM_FILE', help="Full path to sorted BAM file.")
argParser.add_argument('OUTPREFIX', help="Full path to output prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-se', '--single_end', dest="SINGLE_END", help="Specifies whether BAM file contains single-end reads.",action='store_true')
args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## GENERATE INSERT SIZE METRICS           ##
############################################
############################################

OUTDIR = dirname(args.OUTPREFIX)
OUTFILE = join('%s.strand.txt' % (args.OUTPREFIX))

if exists(args.BAM_FILE) and not exists(OUTFILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Determining strand specificity.')
    countDict = uF.getStrandednessFromBAM(args.BAM_FILE,args.SINGLE_END)
    
    senseProp = countDict['sense']/float(countDict['antisense']+countDict['sense'])
    antiSenseProp = 1-senseProp
    forwardProb = 0.5
    if antiSenseProp >= 0.8:
        forwardProb = 0
    elif antiSenseProp <= 0.2:
        forwardProb = 1
    
    uF.makedir(OUTDIR)
    fout = open(OUTFILE,'w')    
    fout.write('RSEM_forward_prob\t%s\n\n' % (forwardProb))
    fout.write('sense_fraction\t%s\n' % (senseProp))
    fout.write('antisense_fraction\t%s\n\n' % (antiSenseProp))
    for pType in ['sense','antisense','discordant','total']:
        fout.write('%s\t%s\n' % (pType,countDict[pType]))
    fout.close()
    
    CompleteFile = join(OUTDIR,'complete/','%s.getStrandSpecificFromBAM.complete' % (basename(args.OUTPREFIX)))
    uF.makedir(dirname(CompleteFile))
    system('touch %s' % (CompleteFile))
    
############################################
############################################
############################################
############################################
