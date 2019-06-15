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

Description = 'Align FastQ file(s) with rsem to transcriptome.'
Epilog = """Example usage: python runRSEMStar.py <CONFIG_FILE> <FASTQ_FILE1> <GENOME_INDEX> <OUTPREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('FASTQ_FILE1', help="Full path to read 1 FASTQ file.")
argParser.add_argument('GENOME_INDEX', help="Full path to RSEM Star genome index including file prefix.")
argParser.add_argument('OUTPREFIX', help="Full path to output prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-fq', '--fastq_file2', type=str, dest="FASTQ_FILE2", default='', help="Full path to read 2 FASTQ file.")
argParser.add_argument('-sf', '--strandedness_file', type=str, dest="STRANDEDNESS_FILE", default='', help="Full path to file containing output of 'getStrandSpecificFromBAM.py'. To manually override can provide a file with two tab-delimited columns with the second column being the --forward-prob value")
argParser.add_argument('-md', '--mark_duplicates', dest="MARK_DUPLICATES", help="Run Picard MarkDuplicates on genome BAM file.",action='store_true')
argParser.add_argument('-ds', '--dedup_bam_suffix', type=str, dest="DEDUP_BAM_SUFFIX", default='markDup', help="Suffix for MarkDuplicates BAM file.")
argParser.add_argument('-nt', '--threads', type=int, dest="NUM_THREADS", default=6, help="Number of threads (default: 6).")
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
SAMTOOLS_EXE = configParser.get('SOFTWARE', 'samtools_exe')
JAVA_18_EXE = configParser.get('SOFTWARE', 'java_18_exe')
PICARD_EXE = configParser.get('SOFTWARE', 'picard_exe')

RSEM_PARAMS = configParser.get('PARAMETERS', 'rsem_params')
JAVA_PARAMS = configParser.get('PARAMETERS', 'java_params')
PICARD_PARAMS = configParser.get('PARAMETERS', 'picard_params')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN,RSEM_BIN,STAR_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## PARAMETERS                             ##
############################################
############################################

PARAM_DICT = {'--forward-prob':0.5}
if exists(args.STRANDEDNESS_FILE):
    fin = open(args.STRANDEDNESS_FILE,'r')
    field,val = fin.readline().strip().split('\t')
    PARAM_DICT['--forward-prob'] = val
    fin.close()

## PAIRED-END PARAMETERS
isPairedEnd = False
if exists(args.FASTQ_FILE1) and exists(args.FASTQ_FILE2):
    PARAM_DICT['--paired-end'] = ''
    isPairedEnd = True
    
## ZIPPED FILES
if args.FASTQ_FILE1[-3:] == '.gz':
    PARAM_DICT['--star-gzipped-read-file'] = ''

## ALIGNER PARAMETERS
PARAM_DICT[' --star '] = ''
PARAM_DICT[' --star-path '] = STAR_BIN

for param,val in sorted(PARAM_DICT.items()):
    if RSEM_PARAMS.find(param) == -1:
        RSEM_PARAMS += ' %s %s' % (param.strip(),val)

############################################
############################################
## ALIGN WITH RSEM STAR                   ##
############################################
############################################

OUTDIR = dirname(args.OUTPREFIX)
TMPDIR = join(OUTDIR,basename(args.OUTPREFIX))

TRANSCRIPT_BAM_FILE = join('%s.transcript.bam' % (args.OUTPREFIX))
TRANSCRIPT_SORTED_BAM_FILE = join('%s.transcriptome.sorted.bam' % (args.OUTPREFIX))
TRANSCRIPT_SORTED_BAM_FLAGSTAT_FILE = join('%s.transcriptome.sorted.bam.flagstat' % (args.OUTPREFIX))

GENOME_BAM_FILE = join('%s.genome.bam' % (args.OUTPREFIX))
if RSEM_PARAMS.find('--star-output-genome-bam') != -1:
    GENOME_BAM_FILE = join('%s.STAR.genome.bam' % (args.OUTPREFIX))
RG_GENOME_BAM_FILE = join('%s.genome.RG.bam' % (args.OUTPREFIX))
GENOME_SORTED_BAM_FILE = join('%s.genome.sorted.bam' % (args.OUTPREFIX))
GENOME_SORTED_BAM_FLAGSTAT_FILE = join('%s.genome.sorted.bam.flagstat' % (args.OUTPREFIX))
GENOME_MARKDUP_SORTED_BAM_FILE = join('%s.genome.%s.sorted.bam' % (args.OUTPREFIX,args.DEDUP_BAM_SUFFIX))
GENOME_MARKDUP_SORTED_BAM_FLAGSTAT_FILE = join('%s.genome.%s.sorted.bam.flagstat' % (args.OUTPREFIX,args.DEDUP_BAM_SUFFIX))

PLOT_FILE = join('%s_model.pdf' % (args.OUTPREFIX))
MARKDUP_METRICS_FILE = join(OUTDIR,'picardMetrics/','%s_MarkDuplicates.metrics.txt' % (basename(args.OUTPREFIX)))

RSEM_SysOutFile = join(OUTDIR,'sysout/','%s.rsem.sysout' % (basename(args.OUTPREFIX)))
ReadGroup_SysOutFile = join(OUTDIR,'sysout/','%s.AddOrReplaceReadGroups.sysout' % (basename(args.OUTPREFIX)))
MarkDup_SysOutFile = join(OUTDIR,'sysout/','%s.MarkDuplicates.sysout' % (basename(args.OUTPREFIX)))

## ALIGN
CommandList = []
if exists(args.FASTQ_FILE1) and not exists(GENOME_SORTED_BAM_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running RSEM.')
    Command = tW.rsemCalculateExpression(RSEM_BIN=RSEM_BIN,FastQFiles1=[args.FASTQ_FILE1],RSEMTranscriptIndex=args.GENOME_INDEX,OutPrefix=args.OUTPREFIX,SysOutFile=RSEM_SysOutFile,FastQFiles2=[args.FASTQ_FILE2],NumThreads=args.NUM_THREADS,RSEMCalcExprParams=RSEM_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Writing RSEM model to file.')
    Command = tW.rsemPlotModel(RSEM_BIN=RSEM_BIN,RSEMFileOutPrefix=args.OUTPREFIX,PDFPlotFile=PLOT_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

## ADD READ GROUPS
if exists(GENOME_BAM_FILE) and not exists(RG_GENOME_BAM_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Adding read groups to Genome BAM.')
    Command = tW.AddOrReplaceReadGroups(JAVA_EXE=JAVA_18_EXE,JAVA_PARAMS=JAVA_PARAMS,PICARD_EXE=PICARD_EXE,BAMFileIn=GENOME_BAM_FILE,BAMFileOut=RG_GENOME_BAM_FILE,RGSampleName=basename(args.OUTPREFIX),SysOutFile=ReadGroup_SysOutFile,RGID=1,Library=1,Platform='illumina',PlatformUnit=1,SeqCentre='null',Description='null',RunDate='null',SORT_ORDER='null',PicardParams=PICARD_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

## SORT GENOME BAM
NUM_READS_IN_FASTQ = -1
if exists(GENOME_BAM_FILE) and not exists(GENOME_SORTED_BAM_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Sorting Genome BAM by coordinate.')
    Command = tW.SortBAMByCoordinate(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFileIn=RG_GENOME_BAM_FILE,BAMFileOut=GENOME_SORTED_BAM_FILE,NumThreads=args.NUM_THREADS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Indexing Genome BAM.')
    Command = tW.indexBAM(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=GENOME_SORTED_BAM_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running flagstat on Genome BAM.')
    Command = tW.flagstat(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=GENOME_SORTED_BAM_FILE,FlagStatFile=GENOME_SORTED_BAM_FLAGSTAT_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

##    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Counting number of reads in FASTQ file.')
##    NUM_READS_IN_FASTQ = uF.numLinesInFile(args.FASTQ_FILE1)
    
##    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Verifying reads in Genome BAM with flagstat output.')
##    isValidBAM = uF.flagStatSTARGenomeBAMValidate(FlagStatFile=GENOME_SORTED_BAM_FLAGSTAT_FILE,NumReadsInFastQ=NUM_READS_IN_FASTQ,isPairedEnd=isPairedEnd)
##    Command = 'touch %s' % (join('%s.fail' % (GENOME_SORTED_BAM_FLAGSTAT_FILE)))
##    if isValidBAM:
##        Command = 'touch %s' % (join('%s.pass' % (GENOME_SORTED_BAM_FLAGSTAT_FILE)))
##    system(Command)
##    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
    
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Deleting unsorted and read group Genome BAM.')
    Command = 'rm %s %s' % (GENOME_BAM_FILE,RG_GENOME_BAM_FILE)
    system(Command)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

## SORT TRANSCRIPTOME BAM
if exists(TRANSCRIPT_BAM_FILE) and not exists(TRANSCRIPT_SORTED_BAM_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Sorting Transcript BAM by coordinate.')
    Command = tW.SortBAMByCoordinate(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFileIn=TRANSCRIPT_BAM_FILE,BAMFileOut=TRANSCRIPT_SORTED_BAM_FILE,NumThreads=args.NUM_THREADS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Indexing Transcript BAM.')
    Command = tW.indexBAM(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=TRANSCRIPT_SORTED_BAM_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running flagstat on Transcript BAM.')
    Command = tW.flagstat(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=TRANSCRIPT_SORTED_BAM_FILE,FlagStatFile=TRANSCRIPT_SORTED_BAM_FLAGSTAT_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
    
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Deleting unsorted Transcript BAM.')
    Command = 'rm %s' % (TRANSCRIPT_BAM_FILE)
    system(Command)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

## MARK DUPLICATES
if args.MARK_DUPLICATES and exists(GENOME_SORTED_BAM_FILE) and not exists(GENOME_MARKDUP_SORTED_BAM_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running Picard MarkDuplicates on Genome BAM.')
    Command = tW.MarkDuplicates(JAVA_EXE=JAVA_18_EXE,JAVA_PARAMS=JAVA_PARAMS,PICARD_EXE=PICARD_EXE,BAMFile=GENOME_SORTED_BAM_FILE,DeDupBAMFile=GENOME_MARKDUP_SORTED_BAM_FILE,MetricsFile=MARKDUP_METRICS_FILE,SysOutFile=MarkDup_SysOutFile,rmDups=False,PicardParams=PICARD_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
                
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Indexing MarkDuplicates Genome BAM.')
    Command = tW.indexBAM(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=GENOME_MARKDUP_SORTED_BAM_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
    
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running flagstat on MarkDuplicates Genome BAM.')
    Command = tW.flagstat(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=GENOME_MARKDUP_SORTED_BAM_FILE,FlagStatFile=GENOME_MARKDUP_SORTED_BAM_FLAGSTAT_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

##    if NUM_READS_IN_FASTQ != -1:
##        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Counting number of reads in FASTQ file.')
##        NUM_READS_IN_FASTQ = uF.numLinesInFile(args.FASTQ_FILE1)
##
##    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Verifying reads in MarkDuplicates BAM with flagstat output.')
##    isValidBAM = uF.flagStatSTARGenomeBAMValidate(FlagStatFile=GENOME_MARKDUP_SORTED_BAM_FLAGSTAT_FILE,NumReadsInFastQ=NUM_READS_IN_FASTQ,isPairedEnd=isPairedEnd)
##    Command = 'touch %s' % (join('%s.fail' % (GENOME_MARKDUP_SORTED_BAM_FLAGSTAT_FILE)))
##    if isValidBAM:
##        Command = 'touch %s' % (join('%s.pass' % (GENOME_MARKDUP_SORTED_BAM_FLAGSTAT_FILE)))
##    system(Command)
##    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

############################################        
############################################
## WRITE COMMAND FILE                     ##
############################################
############################################

if len(CommandList) != 0:
    CompleteFile = join(OUTDIR,'complete/','%s.runRSEMStar.complete' % (basename(args.OUTPREFIX)))
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

##############################################
##############################################
##############################################
##############################################
