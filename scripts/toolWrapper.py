
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

import utilityFunctions as uF

####################################################
####################################################
## RSEM
####################################################
####################################################

def rsemPrepareReference(RSEM_BIN,GenomeFasta,GTFFile,RSEMTranscriptIndex,NumThreads=4,RSEMPrepRefParams=''):

    uF.makedir(dirname(RSEMTranscriptIndex))
    Command = """%srsem-prepare-reference %s --gtf %s --num-threads %s %s %s """ % (RSEM_BIN,RSEMPrepRefParams,GTFFile,NumThreads,GenomeFasta,RSEMTranscriptIndex)
    print Command
    system(Command)
    
    return Command

############################################

def rsemGenerateNGVector(RSEM_BIN,RSEMTranscriptIndex):
    
    Command = """%srsem-generate-ngvector %s.transcripts.fa %s""" % (RSEM_BIN,RSEMTranscriptIndex,RSEMTranscriptIndex)
    print Command
    system(Command)
    
    return Command

############################################

def rsemPlotModel(RSEM_BIN,RSEMFileOutPrefix,PDFPlotFile):

    uF.makedir(dirname(PDFPlotFile))
    Command = """%srsem-plot-model %s %s""" % (RSEM_BIN,RSEMFileOutPrefix,PDFPlotFile)
    print Command
    system(Command)
    
    return Command

############################################

def rsemExtractTranscriptsFromGTF(RSEM_BIN,GTFFile,GenomeFasta,OutPrefix):

    uF.makedir(dirname(OutPrefix))
    Command = """%srsem-extract-reference-transcripts %s 0 %s None 0 %s""" % (RSEM_BIN,OutPrefix,GTFFile,GenomeFasta)
    print Command
    system(Command)
    
    return Command

############################################

def rsemCalculateExpression(RSEM_BIN,FastQFiles1,RSEMTranscriptIndex,OutPrefix,SysOutFile,FastQFiles2=[],NumThreads=4,RSEMCalcExprParams=''):
    
    OutDir = dirname(OutPrefix)
    tTMPDir = join(OutDir,'rsem_tmp/')
    TMPDir = join(tTMPDir,'%s/' % (basename(OutPrefix)))
    uF.makedir(tTMPDir)
    uF.makedir(dirname(SysOutFile))
    Command = """%srsem-calculate-expression --num-threads %s --temporary-folder %s %s %s %s %s %s > %s 2>&1""" % (RSEM_BIN,NumThreads,TMPDir,RSEMCalcExprParams,','.join([x for x in FastQFiles1 if exists(x)]),','.join([x for x in FastQFiles2 if exists(x)]),RSEMTranscriptIndex,OutPrefix,SysOutFile)                                                                           
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))
    
    return Command

####################################################
####################################################
## STAR
####################################################
####################################################

## INPUT  = GenomeFasta
## OUTPUT = OutDir

def StarGenomeIndex(STAR_BIN,GenomeFasta,OutDir,NumThreads=4,StarGenomeIndexParams=''):

    uF.makedir(OutDir)
    Command = """%sSTAR --runMode genomeGenerate --runThreadN %s --genomeFastaFiles %s --genomeDir %s %s""" % (STAR_BIN,NumThreads,GenomeFasta,OutDir,StarGenomeIndexParams)
    print Command
    system(Command)
    
    return Command

####################################################

## INPUT  = FastQFiles, GenomeIndexDir
## OUTPUT = OutDir

def STARAlignReads(STAR_BIN,FastQFiles1,GenomeIndexDir,OutPrefix,SysOutFile,FastQFiles2=[],NumThreads=4,STARAlignReadsParams=''):

    OutDir = dirname(OutPrefix)
    tTMPDir = join(OutDir,'star_tmp/')
    TMPDir = join(tTMPDir,'%s/' % (basename(OutPrefix).strip('.')))
    uF.makedir(tTMPDir)
    uF.makedir(dirname(SysOutFile))
    Command = """%sSTAR --runMode alignReads --readFilesIn %s %s --genomeDir %s --outFileNamePrefix %s --outTmpDir %s --runThreadN %s %s >> %s 2>&1""" % (STAR_BIN,','.join([x for x in FastQFiles1 if exists(x)]),','.join([x for x in FastQFiles2 if exists(x)]),GenomeIndexDir,OutPrefix,TMPDir,NumThreads,STARAlignReadsParams,SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))
    
    return Command

####################################################
####################################################
## QC
####################################################
####################################################

## INPUT  = InputFile
## OUTPUT = OutDir

def FastQC(FASTQC_EXE,JAVA_EXE,InputFile,OutDir,FileFormat='fastq',NumThreads=6):

    uF.makedir(OutDir)
    Command = "%s -q --extract --outdir %s -f %s -t %s --java %s %s" % (FASTQC_EXE,OutDir,FileFormat,NumThreads,JAVA_EXE,InputFile)
    print Command
    system(Command)

    return Command

####################################################

def FastQScreen(FASTQ_SCREEN_EXE,InputFile,OutDir,ConfigFile,SysOutFile,Subset=200000,NumThreads=6):

    uF.makedir(OutDir)
    uF.makedir(dirname(SysOutFile))
    Command = "%s --outdir %s --subset %s --conf %s --threads %s --aligner bowtie2 %s > %s 2>&1" % (FASTQ_SCREEN_EXE,OutDir,Subset,ConfigFile,NumThreads,InputFile,SysOutFile)
    print Command
    system(Command)

    return Command

####################################################

def RNASeqC(JAVA17_EXE,JAVA_PARAMS,RNASEQC_EXE,SamplesFile,OutDir,GenomeFasta,GTFFile,SysOutFile,RNASeqCParams=''):
        
    uF.makedir(OutDir)
    uF.makedir(dirname(SysOutFile))
    Command = """%s %s -jar %s -s %s -o %s -r %s -t %s -gatkFlags "-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY" %s > %s 2>&1 """ % (JAVA17_EXE,JAVA_PARAMS,RNASEQC_EXE,SamplesFile,OutDir,GenomeFasta,GTFFile,RNASeqCParams,SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))
    
    return Command

####################################################
####################################################
## TRIMMING
####################################################
####################################################

def TrimGalore(TRIM_GALORE_EXE,CUTADAPT_BIN,OutDir,FastQFile1,SysOutFile,FastQFile2='',ZipOutput=False,TrimGaloreOptions=''):
    
    uF.makedir(OutDir)
    uF.makedir(dirname(SysOutFile))
   
    if exists(FastQFile2):
        if TrimGaloreOptions.find('--paired') == -1:
            TrimGaloreOptions += ' --paired'

    zipParam = ' --dont_gzip '
    if ZipOutput:
        zipParam = ' --gzip '
    if TrimGaloreOptions.find(zipParam) == -1:
        TrimGaloreOptions += zipParam
    
    Command = "%s --path_to_cutadapt %scutadapt --output_dir %s %s %s >> %s 2>&1" % (TRIM_GALORE_EXE,CUTADAPT_BIN,OutDir,TrimGaloreOptions,' '.join([x for x in [FastQFile1,FastQFile2] if exists(x)]),SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))

    return Command

####################################################

def CutAdapt(CUTADAPT_BIN,OutPrefix,FastQFile1,SysOutFile,FastQFile2='',ZipOutput=False,CutAdaptOptions=''):
    
    uF.makedir(dirname(OutPrefix))
    uF.makedir(dirname(SysOutFile))
    FastQOutFile1 = OutPrefix.strip()+'_1.fastq'
    FastQOutFile2 = OutPrefix.strip()+'_2.fastq'
    if ZipOutput:
        FastQOutFile1 += '.gz'
        FastQOutFile2 += '.gz'
    FastQFiles = [FastQFile1]
    if exists(FastQFile2):
        if CutAdaptOptions.find('--paired-output') == -1:
            CutAdaptOptions += ' --paired-output=%s' % (FastQOutFile2)
        FastQFiles.append(FastQFile2)
    else:
        FastQOutFile1 = OutPrefix.strip()+'.fastq'
        if ZipOutput:
            FastQOutFile1 += '.gz'
    Command = '%scutadapt --output=%s %s %s >> %s 2>&1' % (CUTADAPT_BIN,FastQOutFile1,CutAdaptOptions,' '.join([x for x in [FastQFile1,FastQFile2] if exists(x)]),SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))
    

    return Command

####################################################
####################################################
## SAMTOOLS
####################################################
####################################################

def faidx(SAMTOOLS_EXE,GenomeFasta):

    Command = "%s faidx %s" % (SAMTOOLS_EXE,GenomeFasta)
    print Command
    system(Command)

    return Command

####################################################

def indexBAM(SAMTOOLS_EXE,BAMFile):

    Command = "%s index %s" % (SAMTOOLS_EXE,BAMFile)
    print Command
    system(Command)

    return Command

####################################################

def flagstat(SAMTOOLS_EXE,BAMFile,FlagStatFile):

    uF.makedir(dirname(FlagStatFile))
    Command = "%s flagstat %s > %s" % (SAMTOOLS_EXE,BAMFile,FlagStatFile)
    print Command
    system(Command)

    return Command

####################################################

## INPUT  = BAMFileIn
## OUTPUT = BAMFileOut

def SortBAMByCoordinate(SAMTOOLS_EXE,BAMFileIn,BAMFileOut,NumThreads=4):

    if NumThreads > 6:
        NumThreads = 6
    uF.makedir(dirname(BAMFileOut))
    Command = """%s sort -@ %s %s %s""" % (SAMTOOLS_EXE,NumThreads,BAMFileIn,BAMFileOut[:-4])
    print Command
    system(Command)
    
    return Command

####################################################
####################################################
## PICARD
####################################################
####################################################

## INPUT = GenomeFasta
## OUTPUT = DictFile

def CreateSequenceDictionary(JAVA_EXE,JAVA_PARAMS,PICARD_EXE,GenomeFasta,DictFile,PicardParams=''):

    uF.makedir(dirname(DictFile))
    Command = "%s %s -jar %s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s %s" % (JAVA_EXE,JAVA_PARAMS,PICARD_EXE,GenomeFasta,DictFile,PicardParams)
    print Command
    system(Command)

    return Command

####################################################

## SOURCE = https://www.broadinstitute.org/gatk/events/slides/1506/GATKwr8-B-1-Mapping_and_processing.pdf
## INPUT  = BAMFile
## OUTPUT = DeDupBAMFile, MetricsFile

def MarkDuplicates(JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFile,DeDupBAMFile,MetricsFile,SysOutFile,rmDups=False,PicardParams=''):

    uF.makedir(dirname(DeDupBAMFile))
    uF.makedir(dirname(MetricsFile))
    uF.makedir(dirname(SysOutFile))
    Command = "%s %s -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s %s" % (JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFile,DeDupBAMFile,MetricsFile,PicardParams)
    if rmDups:
        Command += " REMOVE_DUPLICATES=true >> %s 2>&1" % (SysOutFile)
    else:
        Command += " REMOVE_DUPLICATES=false >> %s 2>&1" % (SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))

    return Command

####################################################

def AddOrReplaceReadGroups(JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFileIn,BAMFileOut,RGSampleName,SysOutFile,RGID=1,Library=1,Platform='illumina',PlatformUnit=1,SeqCentre='null',Description='null',RunDate='null',SORT_ORDER='null',PicardParams=''):

    uF.makedir(dirname(BAMFileOut))
    uF.makedir(dirname(SysOutFile))
    Command = "%s %s -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s SORT_ORDER=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s RGDS=%s RGDT=%s %s >> %s 2>&1" % (JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFileIn,BAMFileOut,SORT_ORDER,RGID,Library,Platform,PlatformUnit,RGSampleName,SeqCentre,Description,RunDate,PicardParams,SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))
    
    return Command

####################################################

## INPUT  = BAMFile, GenomeFasta
## OUTPUT = MetricsFile

def CollectMultipleMetrics(JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFile,GenomeFasta,MetricsFile,SysOutFile,PicardParams=''):

    uF.makedir(dirname(MetricsFile))
    uF.makedir(dirname(SysOutFile))
    Command = "%s %s -jar %s CollectMultipleMetrics INPUT=%s OUTPUT=%s REFERENCE_SEQUENCE=%s %s >> %s 2>&1" % (JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFile,MetricsFile,GenomeFasta,PicardParams,SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))

    return Command

####################################################

## INPUT  = BAMFileIn, 
## OUTPUT = BAMFileOut, HistogramFile

def CollectInsertSizeMetrics(JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFileIn,MetricsFileOut,HistogramFile,SysOutFile,PicardParams=''):

    uF.makedir(dirname(MetricsFileOut))
    uF.makedir(dirname(HistogramFile))
    uF.makedir(dirname(SysOutFile))
    Command = "%s %s -jar %s CollectInsertSizeMetrics INPUT=%s OUTPUT=%s HISTOGRAM_FILE=%s %s >> %s 2>&1" % (JAVA_EXE,JAVA_PARAMS,PICARD_EXE,BAMFileIn,MetricsFileOut,HistogramFile,PicardParams,SysOutFile)
    print Command
    system(Command)
    if exists(SysOutFile):
        with open(SysOutFile,'a') as f: f.write('\n\n%s\n\n' % (Command))

    return Command

####################################################
####################################################
####################################################
####################################################

##    ## RSEM COUNTS MATRIX
##    CountsMatrix = '%s%s.%s.counts.matrix' % (OutDir,OutPrefix,entity)
##    if not os.path.exists(CountsMatrix):
##        os.system('%sutil/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix %s %s' % (PEXE.TRINITY_DIR,OutDir + OutPrefix,' '.join(RSEMFileList)))

