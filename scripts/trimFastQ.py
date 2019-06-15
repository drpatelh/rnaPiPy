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

Description = 'Remove adapter and quality trim with TrimGalore or quality trim only with cutadapt.'
Epilog = """Example usage: python trimFastQ.py <CONFIG_FILE> <FASTQ_FILE1> <OUTPREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Software configuration file.")
argParser.add_argument('FASTQ_FILE1', help="Full path to read 1 FASTQ file.")
argParser.add_argument('OUTPREFIX', help="Full path to output prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-fq', '--fastq_file2', type=str, dest="FASTQ_FILE2", default='', help="Full path to read 2 FASTQ file.")
argParser.add_argument('-zo', '--zip_output', dest="ZIP_OUTPUT", help="Zip output fastq file(s).",action='store_true')
argParser.add_argument('-as', '--skip_adapter_removal', dest="SKIP_ADAPTER_REMOVAL", help="Do not run TrimGalore to remove adapter sequences. Quality trimming will still be performed.",action='store_true')
args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')
TRIM_GALORE_EXE = configParser.get('SOFTWARE', 'trim_galore_exe')
CUTADAPT_BIN = configParser.get('SOFTWARE', 'cutadapt_bin')
CUTADAPT_LIB = configParser.get('SOFTWARE', 'cutadapt_lib')

TRIM_GALORE_OPTIONS = configParser.get('PARAMETERS', 'trim_galore_options')
CUTADAPT_QUALTRIM_OPTIONS = configParser.get('PARAMETERS', 'cutadapt_qualtrim_options')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PYTHONPATH"] = '%s:%s' % (':'.join([CUTADAPT_LIB]),environ["PYTHONPATH"])

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## PREPARE OUTPUT FILE NAMES              ##
############################################
############################################

CommandList = []
OUTDIR = dirname(args.OUTPREFIX)

FASTQ_FILEOUT_1 = join('%s_1.fastq' % (args.OUTPREFIX))
FASTQ_FILEOUT_2 = join('%s_2.fastq' % (args.OUTPREFIX))
TrimGalore_FastQFileOut1 = join(OUTDIR,'%s_val_1.fq' % (splitext(basename(args.FASTQ_FILE1))[0]))
TrimGalore_FastQFileOut2 = join(OUTDIR,'%s_val_2.fq' % (splitext(basename(args.FASTQ_FILE2))[0]))

REPORT_FILEOUT_1 = join(OUTDIR,'%s_1.trimming_report.txt' % (args.OUTPREFIX))
REPORT_FILEOUT_2 = join(OUTDIR,'%s_2.trimming_report.txt' % (args.OUTPREFIX))
TrimGalore_Report1 = join(OUTDIR,'%s_trimming_report.txt' % (basename(args.FASTQ_FILE1)))
TrimGalore_Report2 = join(OUTDIR,'%s_trimming_report.txt' % (basename(args.FASTQ_FILE2)))

if args.FASTQ_FILE1[-3:] == '.gz':
    TrimGalore_FastQFileOut1 = join(OUTDIR,'%s_val_1.fq' % (splitext(splitext(basename(args.FASTQ_FILE1))[0])[0]))
    TrimGalore_FastQFileOut2 = join(OUTDIR,'%s_val_2.fq' % (splitext(splitext(basename(args.FASTQ_FILE2))[0])[0]))

if not exists(args.FASTQ_FILE2):
    FASTQ_FILEOUT_1 = join('%s.fastq' % (args.OUTPREFIX))
    REPORT_FILEOUT_1 = join(OUTDIR,'%s.trimming_report.txt' % (args.OUTPREFIX))
    TrimGalore_FastQFileOut1 = join(OUTDIR,'%s_trimmed.fq' % (splitext(basename(args.FASTQ_FILE1))[0]))
    if args.FASTQ_FILE1[-3:] == '.gz':
        TrimGalore_FastQFileOut1 = join(OUTDIR,'%s_trimmed.fq' % (splitext(splitext(basename(args.FASTQ_FILE1))[0])[0]))

if args.ZIP_OUTPUT:
    TrimGalore_FastQFileOut1 += '.gz'
    TrimGalore_FastQFileOut2 += '.gz'
    FASTQ_FILEOUT_1 += '.gz'
    FASTQ_FILEOUT_2 += '.gz'

############################################
############################################
## ADAPTER REMOVAL                        ##
############################################
############################################

if not args.SKIP_ADAPTER_REMOVAL:
    uF.makedir(OUTDIR)
    if not exists(FASTQ_FILEOUT_1) and exists(args.FASTQ_FILE1):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Adapter trimming with TrimGalore.')
        trimgalore_SysOutFile = join(OUTDIR,'sysout/','%s.trimgalore.sysout' % (basename(args.OUTPREFIX)))
        Command = tW.TrimGalore(TRIM_GALORE_EXE=TRIM_GALORE_EXE,CUTADAPT_BIN=CUTADAPT_BIN,OutDir=OUTDIR,FastQFile1=args.FASTQ_FILE1,SysOutFile=trimgalore_SysOutFile,FastQFile2=args.FASTQ_FILE2,ZipOutput=args.ZIP_OUTPUT,TrimGaloreOptions=TRIM_GALORE_OPTIONS)
        CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

        if exists(TrimGalore_FastQFileOut1):
            Command = 'mv %s %s' % (TrimGalore_FastQFileOut1,FASTQ_FILEOUT_1)
            system(Command)
            CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            Command = 'mv %s %s' % (TrimGalore_Report1,REPORT_FILEOUT_1)
            system(Command)
            CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
        if exists(TrimGalore_FastQFileOut2):
            Command = 'mv %s %s' % (TrimGalore_FastQFileOut2,FASTQ_FILEOUT_2)
            system(Command)        
            CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            Command = 'mv %s %s' % (TrimGalore_Report2,REPORT_FILEOUT_2)
            system(Command)
            CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
            
############################################
############################################
## QUALITY TRIMMING                       ##
############################################
############################################

if args.SKIP_ADAPTER_REMOVAL:
    uF.makedir(OUTDIR)
    if not exists(FASTQ_FILEOUT_1) and exists(args.FASTQ_FILE1):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Quality trimming with CutAdapt.')        
        cutadapt_SysOutFile = join(OUTDIR,'sysout/','%s.cutadapt.sysout' % (basename(args.OUTPREFIX)))
        Command = tW.CutAdapt(CUTADAPT_BIN=CUTADAPT_BIN,OutPrefix=args.OUTPREFIX,FastQFile1=args.FASTQ_FILE1,SysOutFile=cutadapt_SysOutFile,FastQFile2=args.FASTQ_FILE2,ZipOutput=args.ZIP_OUTPUT,CutAdaptOptions=CUTADAPT_QUALTRIM_OPTIONS)
        CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

############################################        
############################################
## WRITE COMMAND FILE                     ##
############################################
############################################

if len(CommandList) != 0:
    CompleteFile = join(OUTDIR,'complete/','%s.trimFastQ.complete' % (basename(args.OUTPREFIX)))
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

##############################################
##############################################
##############################################
##############################################
