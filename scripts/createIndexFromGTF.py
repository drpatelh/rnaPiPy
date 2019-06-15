#!/usr/bin/env python

import argparse
import ConfigParser
from time import gmtime, strftime
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import toolWrapper as tW
import utilityFunctions as uF

############################################
## PARSE ARGUMENTS                        ##
############################################

Description = 'Create STAR genome and RSEM STAR transcriptome indices from a gtf file.'
Epilog = """Example usage: python createSTARIndexFromGTF.py <CONFIG_FILE> <GENOME_FASTA_FILE> <GTF_FILE> <OUTDIR> <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('GENOME_FASTA_FILE', help="Full path to genome FASTA file.")
argParser.add_argument('GTF_FILE', help="Full path to GTF file.")
argParser.add_argument('OUTDIR', help="Full path to output directory.")
argParser.add_argument('PREFIX', help="Output prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-nt', '--threads', type=int, dest="NUM_THREADS", default=6, help="Number of threads (default: 6).")
argParser.add_argument('-ri', '--ribosomal_interval_file', type=str, dest="RIBOSOMAL_INTERVAL_FILE", default='',help="Used to create ribosomal LIST file required by RNASeqC. Must tab-delimited file with format 'chr start end strand name' with no header.")
argParser.add_argument('-so', '--star_sjdboverhang', type=int, dest="STAR_SJDBOVERHANG", default=100, help="Only applies if creating RSEM Star index. According to STAR's manual, its ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as well as the ideal value.")
args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')
RSEM_BIN = configParser.get('SOFTWARE', 'rsem_bin')
STAR_BIN = configParser.get('SOFTWARE', 'star_bin')
JAVA_18_EXE = configParser.get('SOFTWARE', 'java_18_exe')
SAMTOOLS_EXE =  configParser.get('SOFTWARE', 'samtools_exe')
PICARD_EXE = configParser.get('SOFTWARE', 'picard_exe')

JAVA_PARAMS = configParser.get('PARAMETERS', 'java_params')
PICARD_PARAMS = configParser.get('PARAMETERS', 'picard_params')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN,RSEM_BIN,STAR_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## CHECK FOR FAI AND DICT FILES           ##
############################################
############################################

GENOME_FILE_DIR = dirname(args.GENOME_FASTA_FILE)
GENOME_PREFIX = splitext(basename(args.GENOME_FASTA_FILE))[0]

GENOME_FAI_FILE = join('%s.fai' % (args.GENOME_FASTA_FILE))
GENOME_DICT_FILE = join(GENOME_FILE_DIR,'%s.dict' % (GENOME_PREFIX))

CommandList = []
if not exists(GENOME_FAI_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating genome FAI file.')
    Command = tW.faidx(SAMTOOLS_EXE=SAMTOOLS_EXE,GenomeFasta=args.GENOME_FASTA_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

if not exists(GENOME_DICT_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating genome DICT file.')
    Command = tW.CreateSequenceDictionary(JAVA_EXE=JAVA_18_EXE,JAVA_PARAMS=JAVA_PARAMS,PICARD_EXE=PICARD_EXE,GenomeFasta=args.GENOME_FASTA_FILE,DictFile=GENOME_DICT_FILE,PicardParams=PICARD_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

############################################
############################################
## CHECK FOR RIBOSMAL LIST FILE           ##
############################################
############################################

if exists(args.RIBOSOMAL_INTERVAL_FILE):
    RIBOSOMAL_LIST_FILE = join(dirname(args.RIBOSOMAL_INTERVAL_FILE),'%s.list' % (splitext(basename(args.RIBOSOMAL_INTERVAL_FILE))[0]))
    if not exists(RIBOSOMAL_LIST_FILE):    
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating rRNA LIST file.')
        Command = 'cat %s %s > %s' % (GENOME_DICT_FILE,args.RIBOSOMAL_INTERVAL_FILE,RIBOSOMAL_LIST_FILE)
        system(Command)
        CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

############################################
############################################
## RSEM STAR                              ##
############################################
############################################

RSEMStarIndexDir = join(args.OUTDIR,'rsem_star/readLen%s/' % (args.STAR_SJDBOVERHANG+1))
RSEMStarExonGTF = join(RSEMStarIndexDir,'%s.exon.gtf' % (args.PREFIX))
RSEMStarTranscriptIndex = join(RSEMStarIndexDir,'%s' % (args.PREFIX))
RSEMStarTranscriptsFasta = join('%s.transcripts.fa' % (RSEMStarTranscriptIndex))
RSEMStarDictFile = join('%s.transcripts.dict' % (RSEMStarTranscriptIndex))

if not exists(join(RSEMStarIndexDir,'sjdbInfo.txt')):

    uF.makedir(RSEMStarIndexDir)
    chdir(RSEMStarIndexDir)

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Filtering GTF file for exon attributes.')
    Command = """awk '$3 == "exon"' %s > %s""" % (args.GTF_FILE,RSEMStarExonGTF)
    system(Command)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
    
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating RSEM Star Index.')
##    Command = tW.rsemPrepareReference(RSEM_BIN=RSEM_BIN,GenomeFasta=args.GENOME_FASTA_FILE,GTFFile=RSEMStarExonGTF,RSEMTranscriptIndex=RSEMStarTranscriptIndex,NumThreads=args.NUM_THREADS,RSEMPrepRefParams='--star --star-path %s --star-sjdboverhang %s --polyA --polyA-length 125' % (STAR_BIN,args.STAR_SJDBOVERHANG))
    Command = tW.rsemPrepareReference(RSEM_BIN=RSEM_BIN,GenomeFasta=args.GENOME_FASTA_FILE,GTFFile=RSEMStarExonGTF,RSEMTranscriptIndex=RSEMStarTranscriptIndex,NumThreads=args.NUM_THREADS,RSEMPrepRefParams='--star --star-path %s --star-sjdboverhang %s' % (STAR_BIN,args.STAR_SJDBOVERHANG))
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating RSEM Star NGVEC file.')
    Command = tW.rsemGenerateNGVector(RSEM_BIN=RSEM_BIN,RSEMTranscriptIndex=RSEMStarTranscriptIndex)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating RSEM Star FAI file.')
    Command = tW.faidx(SAMTOOLS_EXE=SAMTOOLS_EXE,GenomeFasta=RSEMStarTranscriptsFasta)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating RSEM Star DICT file.')
    Command = tW.CreateSequenceDictionary(JAVA_EXE=JAVA_18_EXE,JAVA_PARAMS=JAVA_PARAMS,PICARD_EXE=PICARD_EXE,GenomeFasta=RSEMStarTranscriptsFasta,DictFile=RSEMStarDictFile,PicardParams=PICARD_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

##############################################        
##############################################
#### WRITE COMMAND FILE                     ##
##############################################
##############################################

if len(CommandList) != 0:
    CompleteFile = join(RSEMStarIndexDir,'%s.%s.createIndexFromGTF.complete' % (args.PREFIX,strftime("%d_%m_%Y", gmtime())))
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

##############################################
##############################################
##############################################
##############################################
