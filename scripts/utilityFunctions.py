
import os
import subprocess
import errno
import random
import itertools
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,splitext,exists

import HTSeq
import pysam

####################################################
####################################################
## GENERAL FUNCTIONS
####################################################
####################################################

def makedir(path):
    
    try:
        makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

####################################################

def submitJob(MSUB_EXE,SHFile,NodeReservation,NumThreads,DependCond='afterok',DependJobs=[]):

    MSubStr = '%s -l ' % (MSUB_EXE)
    if len(DependJobs) == 1:
        MSubStr += 'depend=%s:%s,' % (DependCond,DependJobs[0])
    elif len(DependJobs) > 1:
        MSubStr += 'depend=%s:"%s",' % (DependCond,' '.join(DependJobs))
        
    if NodeReservation in ['babs','g7babs96','g7babs60','g7blade']:
        MSubStr += 'nodes=1:%s:ppn=%s %s -N %s -j oe' % (NodeReservation,NumThreads,SHFile,SHFile[:-3])
    else:
        MSubStr += 'nodes=1:ppn=%s %s -N %s -j oe' % (NumThreads,SHFile,SHFile[:-3])

    chdir(dirname(SHFile))
    proc = subprocess.Popen(MSubStr, shell=True, stdout=subprocess.PIPE)
    
    return proc.stdout.read().strip()

####################################################

def numLinesInFile(File):
    
    cmd = "wc -l %s" % (File)
    if File[-3:] == '.gz':
        cmd = "zcat %s | wc -l" % (File)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)

    return int(result.strip().split()[0])

####################################################

def picardInsertMetricsToDict(InsertMetricsFile):

    MetricsDict = {}
    if exists(InsertMetricsFile):
        fin = open(InsertMetricsFile,'r')
        lines = fin.readlines()        
        for idx in range(len(lines)):
            if lines[idx][:len('MEDIAN_INSERT_SIZE')] == 'MEDIAN_INSERT_SIZE':
                MetricsDict = dict(zip(lines[idx].strip().split('\t'),lines[idx+1].strip().split('\t')))
                fin.close()
                break

    return MetricsDict

####################################################

def FAIFileToList(FAIFile):

    FAIList = []
    if exists(FAIFile):
        fin = open(FAIFile,'r')
        for line in fin.readlines():
            lspl = [x.strip() for x in line.strip().split('\t')]
            FAIList.append((lspl[0],int(lspl[1])))
        fin.close()
    
    return FAIList

####################################################

def filterGTFByChromIDs(GTFFileIn,GTFFileOut,ChromIDs):

    if exists(GTFFileIn):
        fin = open(GTFFileIn,'r')
        fout = open(GTFFileOut,'w')
        while True:
            line = fin.readline()
            if line:
                lspl = line.strip().split('\t')
                if lspl[0] in ChromIDs:
                    fout.write(line)
            else:
                fin.close()
                fout.close()
                break
            
####################################################

def getUniqueValsInFileCol(FileIn,colIdx=0,header=False):
    
    uVals = []
    if exists(FileIn):
        fin = open(FileIn,'r')
        if header:
            fin.readline()
        while True:
            line = fin.readline()
            if line:
                lspl = line.strip().split('\t')
                if not lspl[colIdx] in uVals:
                    uVals.append(lspl[colIdx])
            else:
                fin.close()
                break
    
    return uVals
        
####################################################
        
def createPicardListFile(PicardDictFile,IntervalFile,OutFile):

    fin = open(PicardDictFile,'r')
    dictStr = fin.read()
    fin.close()
    fout = open(OutFile,'w')
    fout.write(dictStr)
    fin = open(RibosomalBED,'r')
    nameDict = {}
    while True:
        line = fin.readline()
        if line:
            chrom,start,end,name,score,strand = line.strip().split('\t')
            if not nameDict.has_key(name):
                nameDict[name] = 1
            else:
                nameDict[name] += 1
            fout.write('%s\t%s\t%s\t%s\t%s\n' % (chrom,start,end,strand,name+'.'+str(nameDict[name])))
        else:
            fin.close()
            fout.close()
            break
        
####################################################
####################################################
## SUBSAMPLE FASTQ FUNCTIONS
####################################################
####################################################

def subSampleFastQPE(ReadFraction,FastQFileIn1,FastQFileIn2,FastQFileOut1,FastQFileOut2,Zip=False):

    in1 = iter(HTSeq.FastqReader(FastQFileIn1))
    in2 = iter(HTSeq.FastqReader(FastQFileIn2))
    out1 = open(FastQFileOut1,"w")
    out2 = open(FastQFileOut2,"w")
    
    for read1, read2 in itertools.izip(in1,in2):
       if random.random() < ReadFraction:
          read1.write_to_fastq_file(out1)
          read2.write_to_fastq_file(out2)
          
    in1.close()
    in2.close()
    out1.close()
    out2.close()
    if Zip:
        system('gzip %s' % (FastQFileOut1))
        system('gzip %s' % (FastQFileOut2))

####################################################
  
def subSampleFastQSE(ReadFraction,FastQFileIn,FastQFileOut,Zip=False):

    in1 = iter(HTSeq.FastqReader(FastQFileIn))
    out1 = open(FastQFileOut,"w")
    
    for read1 in in1:
       if random.random() < ReadFraction:
          read1.write_to_fastq_file(out1)

    in1.close()
    out1.close()
    if Zip:
        system('gzip %s' % (FastQFileOut))
        
####################################################
####################################################
## FLAGSTAT FUNCTIONS
####################################################
####################################################

def getNumReadsFromFlagStat(FlagStatFile):
    
    fin = open(FlagStatFile,'r')
    line = fin.readline()
    lspl = line.strip().split()
    fin.close()
    
    return int(lspl[0]) + int(lspl[2])

####################################################

## NEED TO TEST THIS WITH SINGLE-END (DOESNT WORK WITH STAR ANYWAY)

def flagStatSTARGenomeBAMValidate(FlagStatFile,NumReadsInFastQ,isPairedEnd=False):
    
    rv = False
    if os.path.exists(FlagStatFile):
        fin = open(FlagStatFile,'r')
        lines = fin.readlines()
        fin.close()
        if len(lines) > 1:
            count1 = (int(lines[0].split()[0]) - int(lines[1].split()[0]))/2        ## in total - secondary
            count2 = int(lines[5].split()[0])                                       ## paired in sequencing
            count3 = int(lines[6].split()[0])                                       ## read1
            count4 = int(lines[7].split()[0])                                       ## read2
            if isPairedEnd:
                if sum([count1,count2,count3,count4])/6 == NumReadsInFastQ:
                    rv = True
            else:
                if count1 == NumReadsInFastQ:
                    rv = True 
    return rv

####################################################
####################################################
## PYSAM FUNCTIONS
####################################################
####################################################

def getStrandednessFromBAM(BAMFile,isSingleEnd=False):
    ''' Modified version of script provided by Adam. Function calculates the proportion of specificity of reads aligned to a
    sense reference transcriptome. Function takes 2 arguments:
    1)  alignment - A sam/bam file. Sorting unrequired.
    2)  isSingleEnd - boolean whether bam file contains single-end reads.
    '''
    
    count = {'total': 0, 'sense': 0, 'antisense': 0, 'discordant': 0}
    alignFile = pysam.Samfile(BAMFile,"rb")
    for read in alignFile:
        flag = read.flag
        # Skip unmapped reads
        if flag & 4:                        
            continue
        # Count mapped reads
        else:                               
            count['total'] += 1

        strand = 'discordant'
        if isSingleEnd:                     
            if flag & 16:
                strand = 'antisense'
            else:
                strand = 'sense'
        else:              
            # paired and mapped in proper pair
            if flag & 1 and flag & 2:
                # first in pair
                if flag & 64:
                    # read reverse strand and mate forward strand
                    if flag & 16 and not flag & 32:
                         strand = 'antisense'
                    # read forward strand and mate reverse strand
                    elif not flag & 16 and flag & 32:
                        strand = 'sense'
                    
                # second in pair  
                elif flag & 128:
                    # read reverse strand and mate forward strand
                    if flag & 16 and not flag & 32:
                        strand = 'sense'
                    # read forward strand and mate reverse strand
                    elif not flag & 16 and flag & 32:
                        strand = 'antisense'
        count[strand] += 1
    
    return count

####################################################
####################################################
####################################################
####################################################