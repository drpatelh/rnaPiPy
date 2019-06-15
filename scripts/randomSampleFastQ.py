#!/usr/bin/env python

import argparse
import ConfigParser
from time import gmtime, strftime
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import utilityFunctions as uF


## WHAT IF FILE IS GZIPPED?

############################################
## PARSE ARGUMENTS                        ##
############################################

Description = 'Randomly subsample FastQ file(s).'
Epilog = """Example usage: python randomSampleFastQ.py <FASTQ_FILE1> <OUTPREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('FASTQ_FILE1', help="Full path to read 1 FASTQ file.")
argParser.add_argument('OUTPREFIX', help="Full path to output prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-fq', '--fastq_file2', type=str, dest="FASTQ_FILE2", default='', help="Full path to read 2 FASTQ file.")
argParser.add_argument('-ss', '--sample_size', type=int, dest="SAMPLE_SIZE", default=200000, help="Number of reads sampled from FASTQ files (default = 200000).")
argParser.add_argument('-zo', '--zip_output', dest="ZIP_OUTPUT", help="Zip output fastq file(s).",action='store_true')
args = argParser.parse_args()

############################################
############################################
## PREPARE OUTPUT FILE NAMES              ##
############################################
############################################

OUTDIR = dirname(args.OUTPREFIX)
uF.makedir(OUTDIR)

FASTQ_FILEOUT_1 = join('%s_1.fastq' % (args.OUTPREFIX))
FASTQ_FILEOUT_2 = join('%s_2.fastq' % (args.OUTPREFIX))
if not exists(args.FASTQ_FILE2):
    FASTQ_FILEOUT_1 = join('%s.fastq' % (args.OUTPREFIX))

if args.ZIP_OUTPUT:
    FASTQ_FILEOUT_1 += '.gz'
    FASTQ_FILEOUT_2 += '.gz'

############################################
############################################
## RANDOMLY SAMPLE READS                  ##
############################################
############################################

if not exists(FASTQ_FILEOUT_1):
    
    CompleteFile = join(OUTDIR,'complete/','%s.randomSampleFastQ.complete' % (basename(args.OUTPREFIX)))
    uF.makedir(dirname(CompleteFile))
        
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Counting number of reads in FASTQ file.')
    numFastQLines = uF.numLinesInFile(args.FASTQ_FILE1)
    ReadFraction = args.SAMPLE_SIZE/(float(numFastQLines/float(4)))
    
    if exists(args.FASTQ_FILE1) and exists(args.FASTQ_FILE2):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Randomly sampling ~%s reads from PE FASTQ files.' % (args.SAMPLE_SIZE))
        uF.subSampleFastQPE(ReadFraction=ReadFraction,FastQFileIn1=args.FASTQ_FILE1,FastQFileIn2=args.FASTQ_FILE2,FastQFileOut1=FASTQ_FILEOUT_1,FastQFileOut2=FASTQ_FILEOUT_2,Zip=args.ZIP_OUTPUT)
        system('touch %s' % (CompleteFile))
        
    elif exists(args.FASTQ_FILE1) and not exists(args.FASTQ_FILE2):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Randomly sampling ~%s reads from SE FASTQ file.' % (args.SAMPLE_SIZE))
        uF.subSampleFastQSE(ReadFraction=ReadFraction,FastQFileIn=args.FASTQ_FILE1,FastQFileOut=FASTQ_FILEOUT_1,Zip=args.ZIP_OUTPUT)
        system('touch %s' % (CompleteFile))
    
##############################################
##############################################
##############################################
##############################################
