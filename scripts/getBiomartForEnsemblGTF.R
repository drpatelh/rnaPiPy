library(biomaRt);
library(rtracklayer);


## Command line arguements:
## getBiomartForEnsemblGTF.R  species_name  ensembl_release  /path/to/input/ensembl.GTF  /path/to/output/BioMart.annotation
args         <- commandArgs(trailingOnly=TRUE);

species.name <- args[1];
release      <- args[2];
gtf.file     <- args[3];
biomart.file <- args[4];

### Examples:
#
#species.name <- "hsapiens";
#release      <- "75";
#gtf.file     <- "/farm/scratch/rs-bio-lif/mitter01/Projects/RNAseq_pipeline/GTF/Homo_sapiens.GRCh37.75.gtf.gz";
#biomart.file <- "/farm/scratch/rs-bio-lif/mitter01/Projects/RNAseq_pipeline/biomart/Homo_sapiens.GRCh37.75.biomart.txt";
#
#species.name <- "mmusculus";
#release      <- "67";
#gtf.file     <- "/farm/scratch/rs-bio-lif/mitter01/Projects/RNAseq_pipeline/GTF/Mus_musculus.NCBIM37.67.gtf.gz";
#biomart.file <- "/farm/scratch/rs-bio-lif/mitter01/Projects/RNAseq_pipeline/biomart/Mus_musculus.NCBIM37.67.biomart.txt";


## Load GTF file and generate a set of transcript ids
gtf.dat <- import(gtf.file);
tx.ids  <- sort(unique(gtf.dat$transcript_id));
tx.ids  <- tx.ids[!is.na(tx.ids)];


## Query BioMart and write annotation file
if (!file.exists(biomart.file)) {
	dataset    <- paste(species.name,"_gene_ensembl",sep='');

	## Use latest release, regardless of organism
 	host       <- "www.ensembl.org";
	mart       <- "ensembl";
	sym.attr   <- "external_gene_name";

	## Get host URLs for archived releases here:
	## http://www.ensembl.org/info/website/archives/index.html

	## Use archived human (hg19)
	if (species.name == "hsapiens" & release=="75") {
		host     <- "feb2014.archive.ensembl.org";
		marts    <- listMarts(host=host);
                mart     <- marts$biomart[marts$version %in% paste("Ensembl Genes",release,sep=' ')];
		sym.attr <- "hgnc_symbol";
	}

        ## Use archived mouse (mm9)
        if (species.name == "mmusculus" & release=="67") {
                host     <- "may2012.archive.ensembl.org";
                marts    <- listMarts(host=host);
                mart     <- marts$biomart[marts$version %in% paste("Ensembl Genes",release,sep=' ')];
		sym.attr <- "mgi_symbol";
        }


	## Get annotation.  
	ensembl    <- useMart(mart,dataset=dataset,host=host);
        filters    <- listFilters(ensembl);
	
	## Check out attributes for a complete list of available fields to return
        attributes <- listAttributes(ensembl); 
	selected.attributes <- c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","transcript_biotype",sym.attr,"entrezgene","description");
	selected.attributes <- intersect(selected.attributes,attributes$name);
		
	biomart.dat <- getBM(
		attributes = selected.attributes,
		filters    = "ensembl_transcript_id", 
		values     = tx.ids,
		mart       = ensembl);

	if (any(colnames(biomart.dat) %in% c("hgnc_symbol","mgi_symbol"))) {
		colnames(biomart.dat)[colnames(biomart.dat) %in% c("hgnc_symbol","mgi_symbol")] <- "external_gene_name";
	}

	write.table(biomart.dat,file=biomart.file,col.names=TRUE,row.names=FALSE,sep="\t",quote=F);
}


