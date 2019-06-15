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

Description = 'Run RNASeqC.'
Epilog = """Example usage: python runRNASeqC.py <CONFIG_FILE> <BAM_FILES> <OUTDIR>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('BAM_FILES', help="Comma separated list of BAM files.")
argParser.add_argument('GENOME_FASTA', help="Full path to genome FASTA file.")
argParser.add_argument('GTF_FILE', help="Reference gene model in GTF fomat.")
argParser.add_argument('OUTDIR', help="Full path to output directory.")

## OPTIONAL PARAMETERS
argParser.add_argument('-rb', '--ribosomal_list_file', type=str, dest="RIBOSOMAL_LIST_FILE", default='', help="File containing ribosomal intervals in LIST format.")
argParser.add_argument('-se', '--single_end', dest="SINGLE_END", help="Provides RNASeqC with '-singleEnd' parameter. Do not provide if BAM files are a mixture of paired- and single-end data.",action='store_true')
args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')
JAVA_17_EXE = configParser.get('SOFTWARE', 'java_17_exe')
RNASEQC_EXE = configParser.get('SOFTWARE', 'rnaseqc_exe')

RNASEQC_PARAMS = configParser.get('PARAMETERS', 'rnaseqc_params')
JAVA_PARAMS = configParser.get('PARAMETERS', 'java_params')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## ALIGN WITH TOPHAT2                     ##
############################################
############################################

if args.SINGLE_END:
    if RNASEQC_PARAMS.find('-singleEnd') == -1:
        RNASEQC_PARAMS += ' -singleEnd '

SAMPLES_FILE = join(args.OUTDIR,'samples.txt')
REPORT_HTML_FILE = join(args.OUTDIR,'report.html')
RNASeqC_SysOutFile = join(args.OUTDIR,'rnaseqc.sysout')

CommandList = []
if not exists(SAMPLES_FILE) and not exists(REPORT_HTML_FILE):

    ## CREATE SAMPLES FILE
    uF.makedir(args.OUTDIR)
    fout = open(SAMPLES_FILE,'w')
    fout.write('\t'.join(['Sample ID','Bam File','Notes']) + '\n')
    for BAMFile in [x.strip() for x in args.BAM_FILES.split(',')]:
        if exists(BAMFile):
            fout.write('\t'.join([basename(BAMFile)[:-4], BAMFile,'NONE']) + '\n')
    fout.close()

    if exists(args.RIBOSOMAL_LIST_FILE) and RNASEQC_PARAMS.find('-rRNA') == -1:
       RNASEQC_PARAMS += ' -rRNA %s' % (args.RIBOSOMAL_LIST_FILE)
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running RNASeqC.')
    Command = tW.RNASeqC(JAVA17_EXE=JAVA_17_EXE,JAVA_PARAMS=JAVA_PARAMS,RNASEQC_EXE=RNASEQC_EXE,SamplesFile=SAMPLES_FILE,OutDir=args.OUTDIR,GenomeFasta=args.GENOME_FASTA,GTFFile=args.GTF_FILE,SysOutFile=RNASeqC_SysOutFile,RNASeqCParams=RNASEQC_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

############################################        
############################################
## WRITE COMMAND FILE                     ##
############################################
############################################

if len(CommandList) != 0:
    CompleteFile = join(args.OUTDIR,'complete/','runRNASeqC.complete')
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

##############################################
##############################################
##############################################
##############################################
##/farm/babs/redhat6/software/python-2.7.10/bin/python /farm/home/patel35/PYTHON/rnaPiPy/v1.0.0/scripts/runRNASeqC.py /farm/home/patel35/PYTHON/rnaPiPy/v1.0.0/config/softwareConfig.ini /farm/scratch/rs-bio-lif/patel35/projects/rnapipy/hg19/BAM/full/rsem/DOX_R1/DOX_R1.genome.markDup.sorted.bam /farm/home/patel35/PYTHON/rnaPiPy/genome/homo_sapiens/ensembl/GRCh37/release-75/fa/hs.fa /farm/home/patel35/PYTHON/rnaPiPy/genome/homo_sapiens/ensembl/GRCh37/release-75/gtf/hs.gtf / --ribosomal_list_file /farm/home/patel35/PYTHON/rnaPiPy/genome/homo_sapiens/ensembl/GRCh37/release-75/gtf/hs.rRNA.intervals.list
