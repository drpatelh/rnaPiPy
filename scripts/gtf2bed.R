library(rtracklayer);


## Command line arguements:
## getBiomartForEnsemblGTF.R /path/to/input/ensembl.GTF /path/to/output/
args      <- commandArgs(trailingOnly=TRUE);

gtf.file  <- args[1];
out.dir   <- args[2];


## Define output bed files, create out directory if necessary.
gene.file <- paste(out.dir,sub(".gtf.*$","",sub(".*\\/","",gtf.file)),".gene.bed",sep='');
tx.file   <- paste(out.dir,sub(".gtf.*$","",sub(".*\\/","",gtf.file)),".transcript.bed",sep='');
if (!tail(strsplit(out.dir,"")[[1]],n=1) == "/") { out.dir <- paste(out.dir,"/",sep=''); }
if (!file.exists(out.dir)) { dir.create(out.dir,recursive=TRUE); }


## Load GTF file
gtf.dat <- import(gtf.file);


## Transcripts
if (!file.exists(tx.file)) {
	if (any(gtf.dat$type %in% "transcript")) {
		tx.dat  <- gtf.dat[gtf.dat$type %in% "transcript",];
		tx.dat  <- data.frame(
			chr    = as.character(seqnames(tx.dat)),
			start  = start(tx.dat),
			end    = end(tx.dat),
			name   = tx.dat$transcript_id,
			score  = 0,
			strand = strand(tx.dat),
			stringsAsFactors=F);
		write.table(tx.dat,tx.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=F);
	} else {
		tx.exon <- as.data.frame(gtf.dat[gtf.dat$type %in% "exon",]);
		tx.gr   <- GRanges(
			seqnames = as.character(tx.exon$transcript_id),
			ranges   = IRanges(start=as.numeric(tx.exon$start),end=as.numeric(tx.exon$end)),
			strand   = Rle(strand(tx.exon$strand)),
			chr.name = as.character(tx.exon$seqnames));
		tx.gr2   <- split(tx.gr,seqnames(tx.gr));
		tx.range <- as.data.frame(range(tx.gr2));
		tx.range$chr <- tx.exon$seqnames[match(tx.range$seqnames,tx.exon$transcript_id)];
                tx.dat  <- data.frame(
                        chr    = as.character(tx.range$chr),
                        start  = tx.range$start,
                        end    = tx.range$end,
                        name   = tx.range$seqnames,
                        score  = 0,
                        strand = tx.range$strand,
                        stringsAsFactors=F);
                write.table(tx.dat,tx.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=F);
	}
}


## Genes
if (!file.exists(gene.file)) {
	if (any(gtf.dat$type %in% "gene")) {
		gene.dat  <- gtf.dat[gtf.dat$type %in% "gene",];
		gene.dat  <- data.frame(
			chr    = as.character(seqnames(gene.dat)),
			start  = start(gene.dat),
			end    = end(gene.dat),
			name   = gene.dat$gene_id,
			score  = 0,
			strand = strand(gene.dat),
			stringsAsFactors=F);
		write.table(gene.dat,gene.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=F);
	} else {
                tx.exon <- as.data.frame(gtf.dat[gtf.dat$type %in% "exon",]);
                tx.gr   <- GRanges(
                        seqnames = as.character(tx.exon$gene_id),
                        ranges   = IRanges(start=as.numeric(tx.exon$start),end=as.numeric(tx.exon$end)),
                        strand   = Rle(strand(tx.exon$strand)),
                        chr.name = as.character(tx.exon$seqnames));
                tx.gr2   <- split(tx.gr,seqnames(tx.gr));
                tx.range <- as.data.frame(range(tx.gr2));
                tx.range$chr <- tx.exon$seqnames[match(tx.range$seqnames,tx.exon$gene_id)];
                tx.dat  <- data.frame(
                        chr    = as.character(tx.range$chr),
                        start  = tx.range$start,
                        end    = tx.range$end,
                        name   = tx.range$seqnames,
                        score  = 0,
                        strand = tx.range$strand,
                        stringsAsFactors=F);
                write.table(tx.dat,gene.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=F);		
	}
}

