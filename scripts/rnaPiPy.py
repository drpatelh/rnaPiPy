#!/usr/bin/env python

import sys
import argparse
import ConfigParser
from time import gmtime, strftime
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import utilityFunctions as uF

########################################################################################
########################################################################################
############################# PARSE PARAMETERS & DESIGN ################################
########################################################################################
########################################################################################

Description = """RNASeq analysis pipeline. Steps include:
                    \n\t(1) Sample FASTQ file
                    \n\t(2) Align sampled FASTQ file with STAR
                    \n\t(3) Determine strand-specificity from sampled alignment
                    \n\t(4) Adapter/Quality trimming with TrimGalore/Cutadapt
                    \n\t(5) Align and generate counts with RSEM STAR
                    \n\t(6) Genome alignment QC with RNASeqC
                    \n\t(7) FASTQ file QC with FastQC and FastQ Screen
              """
                    ##\n\t(7) Generation of plots for experiment assessment with DESeq2
                    ##\n\t(8) Differential expression with DESeq2                    

Epilog = """Example usage: python rnaPiPy.py -nt 1 -rv any -dr <CONFIG_FILE> <EXPERIMENT_DESIGN_FILE> <OUTDIR> <GENOME_FASTA> <GENOME_INDEX> <ANNOTATION_GTF_FILE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Software configuration file.")
argParser.add_argument('EXPERIMENT_DESIGN_FILE', help="Experiment design file. Must contain following columns ...")   #####
argParser.add_argument('OUTDIR', help="Output directory.")
argParser.add_argument('GENOME_FASTA', help="Genome FASTA file.")
argParser.add_argument('GENOME_INDEX', help="RSEM Star genome index including file prefix.")
argParser.add_argument('ANNOTATION_GTF_FILE', help="Reference gene model in GTF fomat.")

## OPTIONAL PARAMETERS
argParser.add_argument('-rb', '--ribosomal_list_file', type=str, dest="RIBOSOMAL_LIST_FILE", default='', help="File containing ribosomal intervals in LIST format.")
argParser.add_argument('-as', '--skip_adapter_removal', dest="SKIP_ADAPTER_REMOVAL", help="Do not run TrimGalore to remove adapter sequences.",action='store_true')
argParser.add_argument('-ss', '--sample_size', type=int, dest="SAMPLE_SIZE", default=200000, help="Number of reads sampled from FASTQ files (default = 200000).")
argParser.add_argument('-nt', '--threads', type=int, dest="NUM_THREADS", default=6, help="Number of threads (default: 6).")
argParser.add_argument('-rv', '--reservation', type=str, dest="RESERVATION", default='any', help="Cluster reservation. Runs on any free node by default. Can be any, babs, g7babs96, g7babs60, g7blade.")
argParser.add_argument('-ad', '--run_all_dependency', dest="RUN_ALL_DEPENDENCY", help="Run all scripts as job dependencies independently.",action='store_true')
argParser.add_argument('-dr', '--dryrun', dest="DRY_RUN", help="Run script without executing anything on farm. May be useful to check for errors in config file too.",action='store_true')
argParser.add_argument('-ov', '--overwrite', dest="OVERWRITE", help="Rerun analysis and overwrite existing files.",action='store_true')
args = argParser.parse_args()

############################################
############################################
## CHECK ARGUMENTS                        ##
############################################
############################################

for arg in vars(args):
    if arg in ['CONFIG_FILE','EXPERIMENT_DESIGN_FILE','GENOME_FASTA','ANNOTATION_GTF_FILE']:
        if not exists(getattr(args, arg)):
            print 'File does not exists\n%s\n%s' % (arg, getattr(args, arg))
            sys.exit()
    elif arg in ['GENOME_INDEX']:
        if not exists('%s.seq' % (getattr(args, arg))):
            print 'Index does not exists\n%s\n%s' % (arg, getattr(args, arg))
            sys.exit()

############################################
############################################
## PARSE CONFIGURATION FILE               ##
############################################
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

SCRIPT_DIR = configParser.get('SOFTWARE', 'script_dir')
MSUB_EXE = configParser.get('SOFTWARE', 'msub_exe')
PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## PARSE EXPERIMENTAL DESIGN FILE         ##
############################################
############################################

##sample.id	fastq_file1	fastq_file2	group	contrast	insert_mean	insert_stdev
sampleList = []
sampleDesignDict = {}
fin = open(args.EXPERIMENT_DESIGN_FILE)
header = fin.readline().strip().split('\t')
while True:
    line = fin.readline()
    if line:
        lspl = line.strip().split('\t')
        sampleID = lspl[0]
        sampleDesignDict[sampleID] = dict(zip(header,lspl))
        if not exists(sampleDesignDict[sampleID]['fastq_file1']):
            print 'FastQ file does not exists\n%s: %s' % (sampleID,sampleDesignDict[sampleID]['fastq_file1'])
            sys.exit()
        
        sampleDesignDict[sampleID]['isPaired'] = False
        if exists(sampleDesignDict[sampleID]['fastq_file2']):
            sampleDesignDict[sampleID]['isPaired'] = True
        sampleList.append(sampleID)
    else:
        fin.close()
        break

########################################################################################
########################################################################################
############################# SAMPLE LEVEL ANALYSIS ####################################
########################################################################################
########################################################################################

commandGroupList = ['PREP', 'SAMPLE_FASTQ', 'SAMPLE_ALIGN', 'TRIM_FASTQ', 'FULL_ALIGN', 'RNA_QC', 'FASTQ_QC', 'TIDY']

sampleFileDict = {}
sampleCommandDict = {}
for sampleID in sampleList:
    
    ####################################################################################
    ####################################################################################
    ##  DECLARE FILES & DIRECTORIES THAT WILL BE CREATED                              ##
    ####################################################################################
    ####################################################################################
    
    sampleFileDict[sampleID] = {}
    
    RawFastQDir = join(args.OUTDIR,'FASTQ/raw/',sampleID)
    RawFastQPrefix = join(RawFastQDir,sampleID)
    sampleFileDict[sampleID]['RAW_FASTQ_FILE1'] = '%s_1.fastq' % (RawFastQPrefix)
    sampleFileDict[sampleID]['RAW_FASTQ_FILE2'] = '%s_2.fastq' % (RawFastQPrefix)
    if not sampleDesignDict[sampleID]['isPaired']:
        sampleFileDict[sampleID]['RAW_FASTQ_FILE1'] = '%s.fastq' % (RawFastQPrefix)
    if sampleDesignDict[sampleID]['fastq_file1'][-3:] == '.gz':
        sampleFileDict[sampleID]['RAW_FASTQ_FILE1'] += '.gz'
        sampleFileDict[sampleID]['RAW_FASTQ_FILE2'] += '.gz'

    StrandDir = join(args.OUTDIR,'QC/strand/',sampleID)
    SampledDir = join(StrandDir,'sampled')
    SampledPrefix = join(SampledDir,sampleID)
    sampleFileDict[sampleID]['SAMPLED_FASTQ_FILE1'] = '%s_1.fastq' % (SampledPrefix)
    sampleFileDict[sampleID]['SAMPLED_FASTQ_FILE2'] = '%s_2.fastq' % (SampledPrefix)
    if not sampleDesignDict[sampleID]['isPaired']:
        sampleFileDict[sampleID]['SAMPLED_FASTQ_FILE1'] = '%s.fastq' % (SampledPrefix)
    
    sampleFileDict[sampleID]['SAMPLED_GENOME_BAM'] = '%s.genome.sorted.bam' % (SampledPrefix)
    sampleFileDict[sampleID]['SAMPLED_TRANSCRIPTOME_BAM'] = '%s.transcriptome.sorted.bam' % (SampledPrefix)
    
    StrandSpecPrefix = join(StrandDir,sampleID)
    sampleFileDict[sampleID]['SAMPLED_STRANDSPEC_FILE'] = '%s.strand.txt' % (StrandSpecPrefix)
        
    TrimFastQDir = join(args.OUTDIR,'FASTQ/trim/',sampleID)
    TrimFastQPrefix = join(TrimFastQDir,sampleID)
    sampleFileDict[sampleID]['TRIM_FASTQ_FILE1'] = '%s_1.fastq' % (TrimFastQPrefix)
    sampleFileDict[sampleID]['TRIM_FASTQ_FILE2'] = '%s_2.fastq' % (TrimFastQPrefix)
    if not sampleDesignDict[sampleID]['isPaired']:
        sampleFileDict[sampleID]['TRIM_FASTQ_FILE1'] = '%s.fastq' % (TrimFastQPrefix)
    
    BAM_SUFFIX = 'markDup'
    TrimAlignDir = join(args.OUTDIR,'BAM/',sampleID)
    TrimAlignPrefix = join(TrimAlignDir,sampleID)
    sampleFileDict[sampleID]['TRANSCRIPTOME_BAM_FILE'] = '%s.transcriptome.sorted.bam' % (TrimAlignPrefix)
    sampleFileDict[sampleID]['GENOME_BAM_FILE'] = '%s.genome.sorted.bam' % (TrimAlignPrefix)
    sampleFileDict[sampleID]['GENOME_MARKDUPS_BAM_FILE'] = '%s.genome.%s.sorted.bam' % (TrimAlignPrefix,BAM_SUFFIX)
    sampleFileDict[sampleID]['GENE_COUNTS_FILE'] = '%s.genes.results' % (TrimAlignPrefix)
    sampleFileDict[sampleID]['ISOFORM_COUNTS_FILE'] = '%s.isoforms.results' % (TrimAlignPrefix)    
    
    QCDir = join(args.OUTDIR,'QC/')
    RNASeqCDir = join(QCDir,'rnaseqc/',sampleID)
    FastQCDir = join(args.OUTDIR,'QC/fastqc/',sampleID)
    FastQScreenDir = join(args.OUTDIR,'QC/fastq_screen/',sampleID)
    
    CompleteFile = join(args.OUTDIR,'complete/','%s.rnaPiPy.complete' % (sampleID))
    
    ####################################################################################
    ####################################################################################
    ##  OVERWRITE EXISTING FILES IF REQUIRED                                          ##
    ####################################################################################
    ####################################################################################
    
    if args.OVERWRITE:
        for oDir in [RawFastQDir,StrandDir,TrimFastQDir,TrimAlignDir,RNASeqCDir,FastQCDir,FastQScreenDir]:
            if exists(oDir):
                sampleCommandDict[sampleID]['PREP'] += ['rm -rf %s' % (oDir)]
        if exists(CompleteFile):
            sampleCommandDict[sampleID]['PREP'] += ['rm %s' % (CompleteFile)]
    
    ####################################################################################
    ####################################################################################
    ##  GENERATE COMMANDS FOR PIPELINE                                                ##
    ####################################################################################
    ####################################################################################
        
    if not exists(CompleteFile):
        
        sampleCommandDict[sampleID] = dict([(x,[]) for x in commandGroupList])
        
        ############################################
        ############################################
        ## CREATE SOFT-LINK TO RAW FASTQ FILE(S)  ##
        ############################################
        ############################################
        
        if not exists(sampleFileDict[sampleID]['RAW_FASTQ_FILE1']):
            uF.makedir(RawFastQDir)        
            Command = 'ln -s %s %s' % (sampleDesignDict[sampleID]['fastq_file1'],sampleFileDict[sampleID]['RAW_FASTQ_FILE1'])
            sampleCommandDict[sampleID]['PREP'] += [Command]
        if sampleDesignDict[sampleID]['isPaired'] and not exists(sampleFileDict[sampleID]['RAW_FASTQ_FILE2']):
            Command = 'ln -s %s %s' % (sampleDesignDict[sampleID]['fastq_file2'],sampleFileDict[sampleID]['RAW_FASTQ_FILE2'])
            sampleCommandDict[sampleID]['PREP'] += [Command]
            
        ############################################
        ############################################
        ## GENERATE SAMPLED FASTQ FILE(S)         ##
        ############################################
        ############################################
        
        SampledFastQCommand = '%sbin/python %srandomSampleFastQ.py %s %s --sample_size %s' % (PYTHON_DIR,SCRIPT_DIR,sampleFileDict[sampleID]['RAW_FASTQ_FILE1'],SampledPrefix,args.SAMPLE_SIZE)
        if sampleDesignDict[sampleID]['isPaired']:
            SampledFastQCommand += ' --fastq_file2 %s' % (sampleFileDict[sampleID]['RAW_FASTQ_FILE2'])     
        sampleCommandDict[sampleID]['SAMPLE_FASTQ'] += [SampledFastQCommand]
        print '\n' + SampledFastQCommand + '\n'
        
        ############################################
        ############################################
        ## ALIGN SAMPLE FASTQ FILE(S) WITH STAR   ##
        ############################################
        ############################################
    
        SampledAlignCommand = '%sbin/python %srunStar.py %s %s %s %s --threads %s' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['SAMPLED_FASTQ_FILE1'],dirname(args.GENOME_INDEX),SampledPrefix,args.NUM_THREADS)
        if sampleDesignDict[sampleID]['isPaired']:
            SampledAlignCommand += ' --fastq_file2 %s' % (sampleFileDict[sampleID]['SAMPLED_FASTQ_FILE2'])
        sampleCommandDict[sampleID]['SAMPLE_ALIGN'] += [SampledAlignCommand]
        print '\n' + SampledAlignCommand + '\n'
        
        ############################################
        ############################################
        ## DETERMINE STRANDEDNESS                 ##
        ############################################
        ############################################
    
        SampledStrandSpecCommand = '%sbin/python %sgetStrandSpecificFromBAM.py %s %s %s' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['SAMPLED_TRANSCRIPTOME_BAM'],StrandSpecPrefix)
        if not sampleDesignDict[sampleID]['isPaired']:
            SampledStrandSpecCommand += ' --single_end'
        sampleCommandDict[sampleID]['SAMPLE_ALIGN'] += [SampledStrandSpecCommand]
        print '\n' + SampledStrandSpecCommand + '\n'
        
        ############################################
        ############################################
        ## ADAPTER/QUALITY TRIM FASTQ FILE        ##
        ############################################
        ############################################

        trimFastQCommand = '%sbin/python %strimFastQ.py %s %s %s' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['RAW_FASTQ_FILE1'],TrimFastQPrefix)
        if sampleDesignDict[sampleID]['isPaired']:
            trimFastQCommand += ' --fastq_file2 %s ' % (sampleFileDict[sampleID]['RAW_FASTQ_FILE2'])
        if args.SKIP_ADAPTER_REMOVAL:
            trimFastQCommand += ' --skip_adapter_removal'
        sampleCommandDict[sampleID]['TRIM_FASTQ'] += [trimFastQCommand]
        print '\n' + trimFastQCommand + '\n'
        
        ############################################
        ############################################
        ## ALIGN WITH RSEM STAR                   ##
        ############################################
        ############################################
    
        TrimAlignCommand = '%sbin/python %srunRSEMStar.py %s %s %s %s --strandedness_file %s --mark_duplicates --dedup_bam_suffix markDup --threads %s' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['TRIM_FASTQ_FILE1'],args.GENOME_INDEX,TrimAlignPrefix,sampleFileDict[sampleID]['SAMPLED_STRANDSPEC_FILE'],args.NUM_THREADS)
        if sampleDesignDict[sampleID]['isPaired']:
            TrimAlignCommand += ' --fastq_file2 %s' % (sampleFileDict[sampleID]['TRIM_FASTQ_FILE2'])
        sampleCommandDict[sampleID]['FULL_ALIGN'] += [TrimAlignCommand]
        print '\n' + TrimAlignCommand + '\n'

        ############################################
        ############################################
        ## RNASEQC                                ##
        ############################################
        ############################################
    
        RNASeqCCommand = '%sbin/python %srunRNASeqC.py %s %s %s %s %s --ribosomal_list_file %s' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['GENOME_MARKDUPS_BAM_FILE'],args.GENOME_FASTA,args.ANNOTATION_GTF_FILE,RNASeqCDir,args.RIBOSOMAL_LIST_FILE)
        if not sampleDesignDict[sampleID]['isPaired']:
            RNASeqCCommand += ' --single_end'
        sampleCommandDict[sampleID]['RNA_QC']  += [RNASeqCCommand]
        print '\n' + RNASeqCCommand + '\n'
    
        ############################################
        ############################################
        ## FASTQ QC                               ##
        ############################################
        ############################################
    
        if sampleDesignDict[sampleID]['isPaired']:
            FastQCCommand1 = '%sbin/python %sfastQ_QC.py %s %s %s %s --threads %s --skip_fastq_screen' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['TRIM_FASTQ_FILE1'],QCDir,sampleID+'_1',args.NUM_THREADS)
            FastQCCommand2 = '%sbin/python %sfastQ_QC.py %s %s %s %s --threads %s --skip_fastq_screen' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['TRIM_FASTQ_FILE2'],QCDir,sampleID+'_2',args.NUM_THREADS)
            FastQScreenCommand = '%sbin/python %sfastQ_QC.py %s %s %s %s --threads %s --sample_size %s --skip_fastqc' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['RAW_FASTQ_FILE1'],QCDir,sampleID,args.NUM_THREADS,args.SAMPLE_SIZE)
            sampleCommandDict[sampleID]['FASTQ_QC'] += [FastQCCommand1]
            sampleCommandDict[sampleID]['FASTQ_QC'] += [FastQCCommand2]
            sampleCommandDict[sampleID]['FASTQ_QC'] += [FastQScreenCommand]
            print '\n' + FastQCCommand1 + '\n'
            print '\n' + FastQCCommand2 + '\n'
            print '\n' + FastQScreenCommand + '\n'
        else:
            FastQCCommand1 = '%sbin/python %sfastQ_QC.py %s %s %s %s --threads %s --skip_fastq_screen' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['TRIM_FASTQ_FILE1'],QCDir,sampleID,args.NUM_THREADS)
            FastQScreenCommand = '%sbin/python %sfastQ_QC.py %s %s %s %s --threads %s --sample_size %s --skip_fastqc' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,sampleFileDict[sampleID]['RAW_FASTQ_FILE1'],QCDir,sampleID,args.NUM_THREADS,args.SAMPLE_SIZE)
            sampleCommandDict[sampleID]['FASTQ_QC'] += [FastQCCommand1]
            sampleCommandDict[sampleID]['FASTQ_QC'] += [FastQScreenCommand]
            print '\n' + FastQCCommand1 + '\n'
            print '\n' + FastQScreenCommand + '\n'

        ############################################
        ############################################
        ## TIDY UP & WRITE COMPLETE FILE          ##
        ############################################
        ############################################
        
        for oDir in [SampledDir]:
            sampleCommandDict[sampleID]['TIDY'] += ["if [ -d %s ]; then\n\trm -rf %s\nfi" % (oDir,oDir)]
            
        for fileType in ['TRIM_FASTQ_FILE1', 'TRIM_FASTQ_FILE2', 'TRANSCRIPTOME_BAM_FILE','GENOME_BAM_FILE']:
            sampleCommandDict[sampleID]['TIDY'] += ["if [ -f %s ]; then\n\trm %s\nfi" % (sampleFileDict[sampleID][fileType],sampleFileDict[sampleID][fileType])]

        sampleCommandDict[sampleID]['TIDY'] += ['touch %s' % (CompleteFile)]
        uF.makedir(dirname(CompleteFile))

########################################################################################
########################################################################################
############################# EXPERIMENT LEVEL ANALYSIS     ############################
########################################################################################
########################################################################################

############################################
############################################
## RUN DESEQ2                             ##
############################################
############################################

        ## PCA, AND OTHER PLOTS
        ## DESEQ2?

########################################################################################
########################################################################################
############################# SETUP JOBS AND RUN #######################################
########################################################################################
########################################################################################

############################################
############################################
## SETUP JOBS BY SAMPLE                   ##
############################################
############################################

commandDependDict = {'PREP':[],
                     'SAMPLE_FASTQ':['PREP'],
                     'SAMPLE_ALIGN':['SAMPLE_FASTQ'],
                     'TRIM_FASTQ':[],
                     'FULL_ALIGN':['TRIM_FASTQ', 'SAMPLE_ALIGN'],
                     'RNA_QC':['FULL_ALIGN'],
                     'FASTQ_QC':['TRIM_FASTQ'],
                     'TIDY':['FASTQ_QC']}

sampleJobIDDict = dict([(x,dict([(x,[]) for x in commandGroupList])) for x in sampleList])

SH_DIR = join(args.OUTDIR,'sh/')
uF.makedir(SH_DIR)
for sampleID in sampleList:
    if args.RUN_ALL_DEPENDENCY:
        commandIdx = 1
        for commandGroup in commandGroupList:
            commandList = sampleCommandDict[sampleID][commandGroup]
            if len(commandList) != 0:
                SHFile = join(SH_DIR,'%s_%s_%s_%s.sh' % (sampleID,strftime("%d_%m_%Y", gmtime()),commandIdx,commandGroup))
                fout = open(SHFile,'w')
                fout.write('\n%s\n\n' % ('\n\n'.join(commandList)))
                fout.close()
                if not args.DRY_RUN:
                    DependJobs = []
                    for cGroup in commandDependDict[commandGroup]:
                        DependJobs += sampleJobIDDict[sampleID][cGroup]
##                    jobid = commandGroup  ##
                    jobid = uF.submitJob(MSUB_EXE=MSUB_EXE,SHFile=SHFile,NodeReservation=args.RESERVATION,NumThreads=args.NUM_THREADS,DependCond='afterok',DependJobs=DependJobs)
                    sampleJobIDDict[sampleID][commandGroup].append(jobid)
            commandIdx += 1
    else:
        commandList = []
        for commandGroup in commandGroupList:
            commandList += sampleCommandDict[sampleID][commandGroup]
##        for command in commandList:     ##
##            system(command)             ##
    
        if len(commandList) != 0:
            SHFile = join(SH_DIR,'%s_%s_rnaPiPy.sh' % (sampleID,strftime("%d_%m_%Y", gmtime())))
            fout = open(SHFile,'w')
            fout.write('\n%s\n\n' % ('\n\n'.join(commandList)))
            fout.close()
            if not args.DRY_RUN:
##                jobid = 'test'  ##
                jobid = uF.submitJob(MSUB_EXE=MSUB_EXE,SHFile=SHFile,NodeReservation=args.RESERVATION,NumThreads=args.NUM_THREADS,DependCond='afterok',DependJobs=[])
                sampleJobIDDict[sampleID][commandGroupList[0]].append(jobid)

## WRITE FILE TO TRACK ALL JOB IDS
if not args.DRY_RUN:
    for sampleID in sampleList:
        allJobIDList = []
        for vList in sampleJobIDDict[sampleID].values():
            for jobid in vList:
                if not jobid in allJobIDList:
                    allJobIDList.append(jobid)
        allJobIDList.sort()
        if len(allJobIDList) != 0:
            AllJobFile = join(SH_DIR,'%s_%s_JOB_ID.txt' % (sampleID,strftime("%d_%m_%Y", gmtime())))
            fout = open(AllJobFile,'w')
            fout.write('%s\n' % ('\n'.join(allJobIDList)))
            fout.close()
            for jobid in allJobIDList:
                print jobid

############################################
############################################
## SETUP JOBS BY EXPERIMENT               ##
############################################
############################################


##############################################
##############################################
##############################################
##############################################
