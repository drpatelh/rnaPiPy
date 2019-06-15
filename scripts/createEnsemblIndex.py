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

Description = 'Download and create RSEM Star indices for a species available from ENSEMBL.'
Epilog = """Example usage: python createIndexForEnsembl.py <CONFIG_FILE> <GENOME_BUILD> <ANNOTATION> <GENOME_FASTA_FILE> <GENOME_BOWTIE2_INDEX> <OUTDIR> <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Full path to configuration file.")
argParser.add_argument('SPECIES_NAME', help="Species name seperated by underscore e.g. saccharomyces_cerevisiae")
argParser.add_argument('NCBI_BUILD', help="NCBI build e.g. ")
argParser.add_argument('ENSEMBL_RELEASE', help="Ensembl release number e.g. 84")
argParser.add_argument('OUTDIR', help="Full path to output directory.")
argParser.add_argument('PREFIX', help="Output prefix.")

## OPTIONAL PARAMETERS
argParser.add_argument('-nt', '--threads', type=int, dest="NUM_THREADS", default=6, help="Number of threads (default: 6).")
argParser.add_argument('-so', '--star_sjdboverhang', type=int, dest="STAR_SJDBOVERHANG", default=100, help="Only applies if creating RSEM Star index. According to STAR's manual, its ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as well as the ideal value.")
argParser.add_argument('-ag', '--append_gtf_file', type=str, dest="APPEND_GTF_FILE", default='', help="Append GTF file to RSEM transcriptome e.g. ERCC spike-ins. GTF file must match FASTA file to be appended. Has to be used together with '--af' option.")
argParser.add_argument('-af', '--append_fasta_file', type=str, dest="APPEND_FASTA_FILE", default='', help="Append FASTA file to RSEM transcriptome e.g. ERCC spike-ins. Sequences must match those in GTF file to be appended. Has to be used together with '--ag' option.")
argParser.add_argument('-fa', '--fasta_ftp_download_link', type=str, dest="FASTA_FTP_DOWNLOAD_LINK", default='NA',help="Manually provide FTP download link for genome. Have to provide this for all fungi and bacterial species except for S. cerevisiae.")
argParser.add_argument('-gt', '--gtf_ftp_download_link', type=str, dest="GTF_FTP_DOWNLOAD_LINK", default='NA',help="Manually provide FTP download link for genome GTF. Have to provide this for all fungi and bacterial species except for S. cerevisiae.")
argParser.add_argument('-gf', '--gff_ftp_download_link', type=str, dest="GFF_FTP_DOWNLOAD_LINK", default='NA',help="Manually provide FTP download link for genome GTF. Have to provide this for all fungi and bacterial species except for S. cerevisiae.")
args = argParser.parse_args()

############################################
## PARSE CONFIGURATION FILE               ##
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

SCRIPT_DIR = configParser.get('SOFTWARE', 'script_dir')
PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')
SAMTOOLS_EXE = configParser.get('SOFTWARE', 'samtools_exe')
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
## DOWNLOAD & PREPARE FASTA               ##
############################################
############################################

CommandList = []

SPECIES_NAME = args.SPECIES_NAME.lower()
SPECIES_NAME = SPECIES_NAME[0].upper() + SPECIES_NAME[1:]

FASTA_DIR = join(args.OUTDIR,args.NCBI_BUILD,'release-%s/fa/' % (args.ENSEMBL_RELEASE))
FASTA_FILE = join(FASTA_DIR,'%s.fa' % (args.PREFIX))

if not exists(FASTA_FILE):

    uF.makedir(FASTA_DIR)
    chdir(FASTA_DIR)
    
    ENSEMBL_FASTA_FTP = 'ftp://ftp.ensembl.org/pub/release-%s/fasta/%s/dna' % (args.ENSEMBL_RELEASE,SPECIES_NAME.lower())
    ENSEMBL_FASTA = '%s.%s.dna.toplevel.fa' % (SPECIES_NAME,args.NCBI_BUILD)
    if args.NCBI_BUILD in ['GRCh38','GRCm38']:
        ENSEMBL_FASTA = '%s.%s.dna.primary_assembly.fa' % (SPECIES_NAME,args.NCBI_BUILD)
    elif args.NCBI_BUILD in ['GRCh37']:
        ENSEMBL_FASTA = '%s.%s.%s.dna.primary_assembly.fa' % (SPECIES_NAME,args.NCBI_BUILD,args.ENSEMBL_RELEASE)
    elif args.NCBI_BUILD in ['NCBIM37']:
        ENSEMBL_FASTA = '%s.%s.%s.dna.toplevel.fa' % (SPECIES_NAME,args.NCBI_BUILD,args.ENSEMBL_RELEASE)
    if args.FASTA_FTP_DOWNLOAD_LINK != 'NA':
        ENSEMBL_FASTA_FTP = dirname(args.FASTA_FTP_DOWNLOAD_LINK)
        ENSEMBL_FASTA = basename(args.FASTA_FTP_DOWNLOAD_LINK)[:-3]
        
    if not exists(ENSEMBL_FASTA):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Downloading Genome FASTA from Ensembl.')
        Command = 'wget %s/%s.gz' % (ENSEMBL_FASTA_FTP,ENSEMBL_FASTA)
        system(Command)
        CommandList.append(Command)
    
        Command = 'gunzip *.gz'
        system(Command)
        CommandList.append(Command)
    
    if exists(args.APPEND_FASTA_FILE):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Appending FASTA to reference FASTA.')
        Command = 'cat %s %s > %s' % (ENSEMBL_FASTA,args.APPEND_FASTA_FILE,FASTA_FILE)
        system(Command)
        CommandList.append(Command)
        
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating FAI for reference.')
        Command = tW.faidx(SAMTOOLS_EXE=SAMTOOLS_EXE,GenomeFasta=FASTA_FILE)
        CommandList.append(Command)

    else:
        Command = 'ln -s %s %s' % (ENSEMBL_FASTA,FASTA_FILE)
        system(Command)
        CommandList.append(Command)
    
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Creating FAI for reference.')
        Command = tW.faidx(SAMTOOLS_EXE=SAMTOOLS_EXE,GenomeFasta=FASTA_FILE)
        CommandList.append(Command)

############################################
############################################
## DOWNLOAD & PREPARE GTF                 ##
############################################
############################################

GTF_DIR = join(args.OUTDIR,args.NCBI_BUILD,'release-%s/gtf/' % (args.ENSEMBL_RELEASE))
GTF_FILE = join(GTF_DIR,'%s.gtf' % (args.PREFIX))
FILT_GTF_FILE = join(GTF_DIR,'%s.filt.gtf' % (args.PREFIX))

if not exists(GTF_FILE):

    uF.makedir(GTF_DIR)
    chdir(GTF_DIR)
    
    ENSEMBL_GTF_FTP = 'ftp://ftp.ensembl.org/pub/release-%s/gtf/%s' % (args.ENSEMBL_RELEASE,SPECIES_NAME.lower())
    ENSEMBL_GTF = '%s.%s.%s.gtf' % (SPECIES_NAME,args.NCBI_BUILD,args.ENSEMBL_RELEASE)
    if args.GTF_FTP_DOWNLOAD_LINK != 'NA':
        ENSEMBL_GTF_FTP = dirname(args.GTF_FTP_DOWNLOAD_LINK)
        ENSEMBL_GTF = basename(args.GTF_FTP_DOWNLOAD_LINK)[:-3]
    
    ENSEMBL_GFF_FTP = 'ftp://ftp.ensembl.org/pub/release-%s/gff3/%s' % (args.ENSEMBL_RELEASE,SPECIES_NAME.lower())
    ENSEMBL_GFF = '%s.%s.%s.gff3' % (SPECIES_NAME,args.NCBI_BUILD,args.ENSEMBL_RELEASE)
    if args.GFF_FTP_DOWNLOAD_LINK != 'NA':
        ENSEMBL_GFF_FTP = dirname(args.GFF_FTP_DOWNLOAD_LINK)
        ENSEMBL_GFF = basename(args.GFF_FTP_DOWNLOAD_LINK)[:-3]
    
    if not exists(ENSEMBL_GTF):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Downloading GTF from Ensembl.')
        Command = 'wget %s/%s.gz' % (ENSEMBL_GTF_FTP,ENSEMBL_GTF)
        system(Command)
        CommandList.append(Command)
    
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Downloading GFF3 from Ensembl.')
        Command = 'wget %s/%s.gz' % (ENSEMBL_GFF_FTP,ENSEMBL_GFF)
        system(Command)
        CommandList.append(Command)
        
        Command = 'gunzip *.gz'
        system(Command)
        CommandList.append(Command)

    if not exists(FILT_GTF_FILE):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Filtering GTF for chromosomes in reference.')
        Command = '%sbin/python %sfilterGTFByFAI.py %s.fai %s %s' % (PYTHON_DIR,SCRIPT_DIR,FASTA_FILE,ENSEMBL_GTF,FILT_GTF_FILE)
        system(Command)
        CommandList.append(Command)
    
    if exists(args.APPEND_GTF_FILE):
        print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Appending GTF to reference GTF.')
        Command = 'cat %s %s > %s' % (FILT_GTF_FILE,args.APPEND_GTF_FILE,GTF_FILE)
        system(Command)
        CommandList.append(Command)
    else:
        Command = 'ln -s %s %s' % (basename(FILT_GTF_FILE),GTF_FILE)
        system(Command)
        CommandList.append(Command)


############################################
############################################
## FILTER FOR TRANSCRIPTS ONLY (RNASEQC)  ##
############################################
############################################

TRANSCRIPT_GTF_FILE = join(GTF_DIR,'%s.txid.gtf' % (args.PREFIX))
if not exists(TRANSCRIPT_GTF_FILE) and exists(GTF_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Filtering GTF for rows with transcript_id value.')
##    Command = """cat %s | awk '$3 == "transcript"' > %s""" % (GTF_FILE,TRANSCRIPT_GTF_FILE)
    Command = """grep -P "^#|transcript_id" %s > %s"""  % (GTF_FILE,TRANSCRIPT_GTF_FILE)
    system(Command)
    CommandList.append(Command)
    
############################################
############################################
## GTF TO BED                             ##
############################################
############################################

BED_DIR = join(args.OUTDIR,args.NCBI_BUILD,'release-%s/bed/' % (args.ENSEMBL_RELEASE))
BED_FILE = join(BED_DIR,'%s.gene.bed' % (args.PREFIX))        
if not exists(BED_FILE) and exists(GTF_FILE):
    uF.makedir(BED_DIR)
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Converting GTF to BED.')
    Command = '%sRscript %sgtf2bed.R %s %s' % (R_BIN,SCRIPT_DIR,GTF_FILE,BED_DIR)
    system(Command)
    CommandList.append(Command)

############################################
############################################
## EXTRACT RIBOSOMAL EXONS FROM GTF       ##
############################################
############################################

RIBOSOMAL_INTERVAL_FILE = join(GTF_DIR,'%s.rRNA.intervals.txt' % (args.PREFIX))
if not exists(RIBOSOMAL_INTERVAL_FILE) and exists(GTF_FILE):
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Extracting Ribosomal Intervals from GTF.')
    Command = """grep 'gene_biotype "rRNA"' %s | awk '$3 == "exon"' | cut -f1,4,5,7,9 | perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on $."; print join "\t", (@F[0,1,2,3], $1)' | sort -k1V -k2n -k3n > %s""" % (GTF_FILE,RIBOSOMAL_INTERVAL_FILE)
    system(Command)
    CommandList.append(Command)

############################################
############################################
## CREATE GENE MAPPINGS FILE              ##
############################################
############################################

GENE_MAPPING_DIR = join(args.OUTDIR,args.NCBI_BUILD,'release-%s/gene_mapping/' % (args.ENSEMBL_RELEASE))
GENE_MAPPING_FILE = join(GENE_MAPPING_DIR,'%s.gene_mapping.txt' % (args.PREFIX))
if not exists(GENE_MAPPING_FILE) and exists(GTF_FILE):
    uF.makedir(GENE_MAPPING_DIR)
    print '%s: %s' % (strftime("%d-%m-%Y %H:%M:%S", gmtime()),'Downloading Gene Mappings for GTF via BioMart.')
    Command = '%sRscript %sgetBiomartForEnsemblGTF.R %s%s %s %s %s' % (R_BIN,SCRIPT_DIR,SPECIES_NAME[0].lower(),SPECIES_NAME.split('_')[1].lower(),args.ENSEMBL_RELEASE,GTF_FILE,GENE_MAPPING_FILE)
    system(Command)
    CommandList.append(Command)

############################################
############################################
## CREATE RSEM STAR INDEX                 ##
############################################
############################################

INDEX_DIR = join(args.OUTDIR,args.NCBI_BUILD,'release-%s/index/%s/' % (args.ENSEMBL_RELEASE,args.PREFIX))
Command = '%sbin/python %screateIndexFromGTF.py %s %s %s %s %s --threads %s --ribosomal_interval_file %s --star_sjdboverhang %s' % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,FASTA_FILE,GTF_FILE,INDEX_DIR,args.PREFIX,args.NUM_THREADS,RIBOSOMAL_INTERVAL_FILE,args.STAR_SJDBOVERHANG)
system(Command)
CommandList.append(Command)

##############################################        
##############################################
#### WRITE COMMAND FILE                     ##
##############################################
##############################################

if len(CommandList) != 0:
    CompleteFile = join(args.OUTDIR,args.NCBI_BUILD,'release-%s/complete/%s.readLen%s.%s.createEnsemblIndex.complete' % (args.ENSEMBL_RELEASE,args.PREFIX,args.STAR_SJDBOVERHANG+1,strftime("%d_%m_%Y", gmtime())))
    uF.makedir(dirname(CompleteFile))
    fout = open(CompleteFile,'w')
    fout.write('\n' + '\n\n'.join(CommandList) + '\n')
    fout.close()

############################################
############################################
############################################
############################################
############################################
############################################
