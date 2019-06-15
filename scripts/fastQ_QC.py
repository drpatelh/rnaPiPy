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

Description = 'Run FastQC and FastQ screen on FastQ file(s).'
Epilog = """Example usage: python fastQ_QC.py <CONFIG_FILE> <FASTQ_FILE> <OUTDIR>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('FASTQ_FILE', help="Full path to FASTQ file.")
argParser.add_argument('OUTDIR', help="Full path to output directory.")
argParser.add_argument('SAMPLE_PREFIX', help="Sample prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-nt', '--threads', type=int, dest="NUM_THREADS", default=6, help="Number of threads (default: 6).")
argParser.add_argument('-fs', '--skip_fastq_screen', dest="SKIP_FASTQ_SCREEN", help="Do not run FastQ Screen.",action='store_true')
argParser.add_argument('-qf', '--skip_fastqc', dest="SKIP_FASTQC", help="Do not run FastQC.",action='store_true')
argParser.add_argument('-ss', '--sample_size', type=int, dest="SAMPLE_SIZE", default=200000, help="Number of reads sampled from FASTQ files (default = 200000).")
args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')
JAVA_18_EXE = configParser.get('SOFTWARE', 'java_18_exe')
FASTQC_EXE = configParser.get('SOFTWARE', 'fastqc_exe')
FASTQ_SCREEN_EXE = configParser.get('SOFTWARE', 'fastq_screen_exe')
FASTQ_SCREEN_CONFIG_FILE = configParser.get('SOFTWARE', 'fastq_screen_config_file')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## RANDOMLY SAMPLE READS                  ##
############################################
############################################

FASTQC_DIR = join(args.OUTDIR,'fastqc/',args.SAMPLE_PREFIX)
FASTQC_HTML_FILE = join(FASTQC_DIR,'%s_fastqc.html' % (args.SAMPLE_PREFIX))
FASTQC_ZIP_FILE = join(FASTQC_DIR,'%s_fastqc.zip' % (args.SAMPLE_PREFIX))

FASTQ_SCREEN_DIR = join(args.OUTDIR,'fastq_screen/',args.SAMPLE_PREFIX)
FASTQ_SCREEN_PNG_FILE = join(args.OUTDIR,'fastq_screen/',args.SAMPLE_PREFIX,'%s_screen.png' % (args.SAMPLE_PREFIX))
FASTQ_SCREEN_TXT_FILE = join(args.OUTDIR,'fastq_screen/',args.SAMPLE_PREFIX,'%s_screen.txt' % (args.SAMPLE_PREFIX))
FastQScreen_SysOutFile = join(args.OUTDIR,'fastq_screen/',args.SAMPLE_PREFIX,'sysout/','%s.fastqscreen.sysout' % (args.SAMPLE_PREFIX))

FastQScreenOutPrefix = splitext(basename(args.FASTQ_FILE))[0]
if args.FASTQ_FILE[-3:] == '.gz':
    FastQScreenOutPrefix = splitext(splitext(basename(args.FASTQ_FILE))[0])[0]
FastQScreenPngFile = join(args.OUTDIR,'fastq_screen/',args.SAMPLE_PREFIX,'%s_screen.png' % (FastQScreenOutPrefix))
FastQScreenTxtFile = join(args.OUTDIR,'fastq_screen/',args.SAMPLE_PREFIX,'%s_screen.txt' % (FastQScreenOutPrefix))

if exists(args.FASTQ_FILE):
    CommandList = []
    if not args.SKIP_FASTQC:
        if not exists(FASTQC_HTML_FILE):
            uF.makedir(FASTQC_DIR)
            print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running FastQC.')
            Command = tW.FastQC(FASTQC_EXE=FASTQC_EXE,JAVA_EXE=JAVA_18_EXE,InputFile=args.FASTQ_FILE,OutDir=FASTQC_DIR,FileFormat='fastq',NumThreads=args.NUM_THREADS)
            CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            
            if exists(FASTQC_ZIP_FILE):
                print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Deleting FastQC ZIP file.')
                Command = 'rm %s' % (FASTQC_ZIP_FILE)
                system(Command)
                CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            
            CompleteFile = join(FASTQC_DIR,'complete/','%s.fastq_QC.complete' % (args.SAMPLE_PREFIX))
            uF.makedir(dirname(CompleteFile))
            fout = open(CompleteFile,'w')
            fout.write('\n' + '\n\n'.join(CommandList) + '\n')
            fout.close()
    
    CommandList = []
    if not args.SKIP_FASTQ_SCREEN:
        if not exists(FASTQ_SCREEN_PNG_FILE):
            uF.makedir(FASTQ_SCREEN_DIR)
            print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running FastQ Screen.')
            Command = tW.FastQScreen(FASTQ_SCREEN_EXE=FASTQ_SCREEN_EXE,InputFile=args.FASTQ_FILE,OutDir=FASTQ_SCREEN_DIR,ConfigFile=FASTQ_SCREEN_CONFIG_FILE,SysOutFile=FastQScreen_SysOutFile,Subset=args.SAMPLE_SIZE,NumThreads=args.NUM_THREADS)
            CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            
            if exists(FastQScreenPngFile) and not exists(FASTQ_SCREEN_PNG_FILE):
                Command = 'mv %s %s' % (FastQScreenPngFile,FASTQ_SCREEN_PNG_FILE)
                system(Command)
                CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            if exists(FastQScreenTxtFile) and not exists(FASTQ_SCREEN_TXT_FILE):
                Command = 'mv %s %s' % (FastQScreenTxtFile,FASTQ_SCREEN_TXT_FILE)
                system(Command)        
                CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

            CompleteFile = join(FASTQ_SCREEN_DIR,'complete/','%s.fastq_QC.complete' % (args.SAMPLE_PREFIX))
            uF.makedir(dirname(CompleteFile))
            fout = open(CompleteFile,'w')
            fout.write('\n' + '\n\n'.join(CommandList) + '\n')
            fout.close()

##############################################
##############################################
##############################################
##############################################
