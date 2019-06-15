#!/usr/bin/env python

import math
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

Description = 'Run Picard CollectInsertSizeMetrics and extract mean and sd of insert size.'
Epilog = """Example usage: python getInsertSizeFromBAM.py <CONFIG_FILE> <BAM_FILE> <OUTPREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('BAM_FILE', help="Full path to sorted BAM file.")
argParser.add_argument('OUTPREFIX', help="Full path to output prefix.")

args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')
JAVA_18_EXE = configParser.get('SOFTWARE', 'java_18_exe')
PICARD_EXE = configParser.get('SOFTWARE', 'picard_exe')

JAVA_PARAMS = configParser.get('PARAMETERS', 'java_params')
PICARD_PARAMS = configParser.get('PARAMETERS', 'picard_params')

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
METRICS_FILE = join('%s.insert_size_metrics' % (args.OUTPREFIX))
HISTOGRAM_FILE = join('%s.insert_size_histogram.pdf' % (args.OUTPREFIX))
SysOutFile = join(OUTDIR,'%s.CollectInsertSizeMetrics.sysout' % (basename(args.OUTPREFIX)))

if exists(args.BAM_FILE) and not exists(METRICS_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running Picard CollectInsertSizeMetrics.')
    Command = tW.CollectInsertSizeMetrics(JAVA_EXE=JAVA_18_EXE,JAVA_PARAMS=JAVA_PARAMS,PICARD_EXE=PICARD_EXE,BAMFileIn=args.BAM_FILE,MetricsFileOut=METRICS_FILE,HistogramFile=HISTOGRAM_FILE,SysOutFile=SysOutFile,PicardParams=PICARD_PARAMS)
    CommandList = ['%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command)]

    CompleteFile = join(OUTDIR,'complete/','%s.getInsertSizeFromBAM.complete' % (basename(args.OUTPREFIX)))
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

if exists(METRICS_FILE):
    insertMetricsDict = uF.picardInsertMetricsToDict(InsertMetricsFile=METRICS_FILE)
    insertMean = int(round(float(insertMetricsDict['MEAN_INSERT_SIZE'])))
    insertStdDev = int(math.ceil(float(insertMetricsDict['STANDARD_DEVIATION'])))
    print '%s\t%s' % (insertMean,insertStdDev)

############################################
############################################
############################################
############################################