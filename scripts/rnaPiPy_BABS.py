#!/usr/bin/env python

import sys
import argparse
import ConfigParser
from time import gmtime, strftime
from os import environ,chdir,system,makedirs
from os.path import basename,dirname,join,splitext,exists

########################################################################################
########################################################################################
############################# PARSE PARAMETERS & DESIGN ################################
########################################################################################
########################################################################################

Description = """Custom RNASeq analysis pipeline for BABS. Steps include:
                    \n\t(1) Sample FASTQ file
                    \n\t(2) Align sampled FASTQ file with STAR
                    \n\t(3) Determine strand-specificity from sampled alignment
                    \n\t(4) Adapter/Quality trimming with TrimGalore/Cutadapt
                    \n\t(5) Align and generate counts with RSEM STAR
                    \n\t(6) Genome alignment QC with RNASeqC
                    \n\t(7) FASTQ file QC with FastQC and FastQ Screen
              """

Epilog = """Example usage: python rnaPiPy_BABS.py -nt 1 -rv any <CONFIG_FILE> <EXPERIMENT_DESIGN_FILE> <SPECIES_NAME> <OUTDIR>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('CONFIG_FILE', help="Software configuration file.")
argParser.add_argument('EXPERIMENT_DESIGN_FILE', help="Experiment design file. Must contain following columns ...")   #####
argParser.add_argument('OUTDIR', help="Output directory.")
argParser.add_argument('SPECIES_NAME', help="Full species name separated by '_' e.g. homo_sapiens")

## OPTIONAL PARAMETERS
argParser.add_argument('-er', '--ercc_spike_in', dest="ERCC_SPIKE_IN", help="Include ERCC spike-ins in transcriptome alignment.",action='store_true')
argParser.add_argument('-as', '--skip_adapter_removal', dest="SKIP_ADAPTER_REMOVAL", help="Do not run TrimGalore to remove adapter sequences.",action='store_true')
argParser.add_argument('-nt', '--threads', type=int, dest="NUM_THREADS", default=6, help="Number of threads (default: 6).")
argParser.add_argument('-rv', '--reservation', type=str, dest="RESERVATION", default='any', help="Cluster reservation. Runs on any free node by default. Can be any, babs, g7babs96, g7babs60, g7blade.")
argParser.add_argument('-ad', '--run_all_dependency', dest="RUN_ALL_DEPENDENCY", help="Run all scripts as job dependencies independently.",action='store_true')
argParser.add_argument('-ov', '--overwrite', dest="OVERWRITE", help="Rerun analysis and overwrite existing files.",action='store_true')
args = argParser.parse_args()

############################################
############################################
## PARSE CONFIGURATION FILE               ##
############################################
############################################

configParser = ConfigParser.SafeConfigParser()
configParser.read(args.CONFIG_FILE)

SCRIPT_DIR = configParser.get('SOFTWARE', 'script_dir')
PYTHON_DIR = configParser.get('SOFTWARE', 'python_dir')
R_BIN = configParser.get('SOFTWARE', 'r_bin')

############################################
## SET ENVIRONMENT VARIABLES              ##
############################################

environ["PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'bin/'),R_BIN]),environ["PATH"])
environ["LD_LIBRARY_PATH"] = '%s:%s' % (':'.join([join(PYTHON_DIR,'lib/'),"/farm/babs/redhat6/lib/"]),environ["LD_LIBRARY_PATH"])

############################################
############################################
## GET GENOME FILE PATHS                  ##
############################################
############################################

GENOME_DIR = '/farm/home/patel35/PYTHON/rnaPiPy/genome/'

SPECIES_NAME = args.SPECIES_NAME.lower()
SPECIES_ABBREV = SPECIES_NAME[0]+ SPECIES_NAME.split('_')[1][0]
if args.ERCC_SPIKE_IN:
    SPECIES_ABBREV += '.ercc'
    
SpeciesDict = {'homo_sapiens':{'build':'GRCh38','release':'84'},
               'mus_musculus':{'build':'GRCm38','release':'84'},
               'gallus_gallus':{'build':'Galgal4','release':'84'},
               'xenopus_tropicalis':{'build':'JGI_4.2','release':'84'},
               'danio_rerio':{'build':'GRCz10','release':'84'},
               'caenorhabditis_elegans':{'build':'WBcel235','release':'84'},
               'drosophila_melanogaster':{'build':'BDGP6','release':'84'},
               'saccharomyces_cerevisiae':{'build':'R64-1-1','release':'84'},
               'schizosaccharomyces_pombe':{'build':'ASM294v2.31','release':'31'},
               'escherichia_coli':{'build':'k12_mg1655','release':'29'}}

if SPECIES_NAME not in SpeciesDict.keys():
    print 'Please provide valid species name.\n'
    sys.exit()

RELEASE_DIR = '%s%s/ensembl/%s/release-%s/' % (GENOME_DIR,SPECIES_NAME,SpeciesDict[SPECIES_NAME]['build'],SpeciesDict[SPECIES_NAME]['release'])
GENOME_FASTA = '%sfa/%s.fa' % (RELEASE_DIR,SPECIES_ABBREV)
GENOME_INDEX = '%sindex/%s/rsem_star/readLen101/%s' % (RELEASE_DIR,SPECIES_ABBREV,SPECIES_ABBREV)
ANNOTATION_GTF_FILE = '%sgtf/%s.txid.gtf' % (RELEASE_DIR,SPECIES_ABBREV)

RIBOSOMAL_LIST_FILE = '%sgtf/%s.rRNA.intervals.list' % (RELEASE_DIR,SPECIES_ABBREV)

############################################
############################################
## RUN RNAPIPY                            ##
############################################
############################################

SAMPLE_SIZE = 200000

Command = """%sbin/python %srnaPiPy.py %s %s %s %s %s %s --ribosomal_list_file %s --sample_size %s --threads %s --reservation %s""" % (PYTHON_DIR,SCRIPT_DIR,args.CONFIG_FILE,args.EXPERIMENT_DESIGN_FILE,args.OUTDIR,GENOME_FASTA,GENOME_INDEX,ANNOTATION_GTF_FILE,RIBOSOMAL_LIST_FILE,SAMPLE_SIZE,args.NUM_THREADS,args.RESERVATION)
if args.SKIP_ADAPTER_REMOVAL:
    Command += ' --skip_adapter_removal'
if args.RUN_ALL_DEPENDENCY:
    Command += ' --run_all_dependency'
if args.OVERWRITE:
    Command += ' --overwrite'
print Command
system(Command)

############################################
############################################
############################################
############################################
