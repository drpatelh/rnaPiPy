#!/usr/bin/env python

import argparse
import ConfigParser
from time import gmtime, strftime
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import toolWrapper as tW
import utilityFunctions as uF

############################################
## NOTES                                  ##
############################################

##GeneralParams = ['--quantMode TranscriptomeSAM', '--outSAMtype BAM Unsorted']
##EncodeSTARParams = ['--outFilterType BySJout', '--outFilterMultimapNmax 20', '--outFilterMismatchNmax 999', '--alignSJoverhangMin 8', '--alignSJDBoverhangMin 1', '--alignIntronMin 20', '--alignIntronMax 1000000',  '--alignMatesGapMax 1000000']
##RSEMSTARParams = ['--outSAMunmapped Within', '--outSAMattributes NH HI AS NM MD', '--outFilterMismatchNoverLmax 0.04', '--sjdbScore 1', '--genomeLoad NoSharedMemory', '--outSAMheaderHD \@HD VN:1.4 SO:unsorted']

############################################
## PARSE ARGUMENTS                        ##
############################################

Description = 'Align FastQ file(s) with Star to genome/transcriptome.'
Epilog = """Example usage: python runStar.py <CONFIG_FILE> <FASTQ_FILE1> <GENOME_INDEX_DIR> <OUTPREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('FASTQ_FILE1', help="Full path to read 1 FASTQ file.")
argParser.add_argument('GENOME_INDEX_DIR', help="Full path to directory containing Star genome index.")
argParser.add_argument('OUTPREFIX', help="Full path to output prefix. Trailing '.' character will be automatically appended.")

## OPTIONAL PARAMETERS
argParser.add_argument('-fq', '--fastq_file2', type=str, dest="FASTQ_FILE2", default='', help="Full path to read 2 FASTQ file.")
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
STAR_BIN = configParser.get('SOFTWARE', 'star_bin')
SAMTOOLS_EXE = configParser.get('SOFTWARE', 'samtools_exe')
JAVA_18_EXE = configParser.get('SOFTWARE', 'java_18_exe')
PICARD_EXE = configParser.get('SOFTWARE', 'picard_exe')

STAR_PARAMS = configParser.get('PARAMETERS', 'star_params')
JAVA_PARAMS = configParser.get('PARAMETERS', 'java_params')
PICARD_PARAMS = configParser.get('PARAMETERS', 'picard_params')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN,STAR_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## PARAMETERS                             ##
############################################
############################################

OUTDIR = dirname(args.OUTPREFIX)

isPairedEnd = False
if exists(args.FASTQ_FILE1) and exists(args.FASTQ_FILE2):
    isPairedEnd = True
    
if args.FASTQ_FILE1[-3:] == '.gz':
    STAR_PARAMS += ' --readFilesCommand zcat'
if STAR_PARAMS.find('--outSAMattrRGline') == -1:
    STAR_PARAMS += ' --outSAMattrRGline ID:1 PU:1 SM:%s LB:1 PL:illumina' % (basename(args.OUTPREFIX))
    
############################################
############################################
## ALIGN WITH RSEM STAR                   ##
############################################
############################################

TRANSCRIPT_BAM_FILE = join('%s.Aligned.toTranscriptome.out.bam' % (args.OUTPREFIX))
TRANSCRIPT_SORTED_BAM_FILE = join('%s.transcriptome.sorted.bam' % (args.OUTPREFIX))
TRANSCRIPT_SORTED_BAM_FLAGSTAT_FILE = join('%s.transcriptome.sorted.bam.flagstat' % (args.OUTPREFIX))

GENOME_BAM_FILE = join('%s.Aligned.out.bam' % (args.OUTPREFIX))
GENOME_SORTED_BAM_FILE = join('%s.genome.sorted.bam' % (args.OUTPREFIX))
GENOME_SORTED_BAM_FLAGSTAT_FILE = join('%s.genome.sorted.bam.flagstat' % (args.OUTPREFIX))

STAR_SysOutFile = join(OUTDIR,'sysout/','%s.star.sysout' % (basename(args.OUTPREFIX)))

GENOME_MARKDUP_SORTED_BAM_FILE = join('%s.genome.%s.sorted.bam' % (args.OUTPREFIX,args.DEDUP_BAM_SUFFIX))
GENOME_MARKDUP_SORTED_BAM_FLAGSTAT_FILE = join('%s.genome.%s.sorted.bam.flagstat' % (args.OUTPREFIX,args.DEDUP_BAM_SUFFIX))
MARKDUP_METRICS_FILE = join(OUTDIR,'picardMetrics/','%s_MarkDuplicates.metrics.txt' % (basename(args.OUTPREFIX)))
MarkDup_SysOutFile = join(OUTDIR,'sysout/','%s.MarkDuplicates.sysout' % (basename(args.OUTPREFIX)))

## ALIGN
CommandList = []
if exists(args.FASTQ_FILE1) and not exists(GENOME_SORTED_BAM_FILE): 
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running Star.')
    Command = tW.STARAlignReads(STAR_BIN=STAR_BIN,FastQFiles1=[args.FASTQ_FILE1],GenomeIndexDir=args.GENOME_INDEX_DIR,OutPrefix=args.OUTPREFIX+'.',SysOutFile=STAR_SysOutFile,FastQFiles2=[args.FASTQ_FILE2],NumThreads=args.NUM_THREADS,STARAlignReadsParams=STAR_PARAMS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

## SORT GENOME BAM
NUM_READS_IN_FASTQ = -1
if exists(GENOME_BAM_FILE) and not exists(GENOME_SORTED_BAM_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Sorting Genome BAM by coordinate.')
    Command = tW.SortBAMByCoordinate(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFileIn=GENOME_BAM_FILE,BAMFileOut=GENOME_SORTED_BAM_FILE,NumThreads=args.NUM_THREADS)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Indexing Genome BAM.')
    Command = tW.indexBAM(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=GENOME_SORTED_BAM_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))

    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Running flagstat on Genome BAM.')
    Command = tW.flagstat(SAMTOOLS_EXE=SAMTOOLS_EXE,BAMFile=GENOME_SORTED_BAM_FILE,FlagStatFile=GENOME_SORTED_BAM_FLAGSTAT_FILE)
    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
    
##    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Counting number of reads in FASTQ file.')
##    NUM_READS_IN_FASTQ = uF.numLinesInFile(args.FASTQ_FILE1)
##    
##    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Verifying reads in Genome BAM with flagstat output.')
##    isValidBAM = uF.flagStatSTARGenomeBAMValidate(FlagStatFile=GENOME_SORTED_BAM_FLAGSTAT_FILE,NumReadsInFastQ=NUM_READS_IN_FASTQ,isPairedEnd=isPairedEnd)
##    Command = 'touch %s' % (join('%s.fail' % (GENOME_SORTED_BAM_FLAGSTAT_FILE)))
##    if isValidBAM:
##        Command = 'touch %s' % (join('%s.pass' % (GENOME_SORTED_BAM_FLAGSTAT_FILE)))
##    system(Command)
##    CommandList.append('%s\n%s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),Command))
 
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Deleting unsorted Genome BAM.')
    Command = 'rm %s' % (GENOME_BAM_FILE)
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
    CompleteFile = join(OUTDIR,'complete/','%s.runStar.complete' % (basename(args.OUTPREFIX)))
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

##############################################
##############################################
##############################################
##############################################
