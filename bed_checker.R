
library(data.table)
library(snpStats, quietly=TRUE)

verbose <- TRUE

args.orig <- commandArgs(trailingOnly=TRUE)

if ('--noprompt' %in% args.orig) {
    verbose <- FALSE
    args.orig <- setdiff(args.orig, '--noprompt')
}

if (length(args.orig) != 2) {
    stop("Usage: Rscript bed_checker.R genotype.file bed.file")
}

genotype.file <- args.orig[1]
bed.file <- args.orig[2]

bed.file <- gsub('.bed$', '', bed.file)
bed.file <- paste0(bed.file, '.bed')
bim.file <- gsub('.bed$', '.bim', bed.file)

cat("\n\n")
cat(paste0("Genotype file is :", genotype.file, "\n"))
cat(paste0("Bed file is: ", bed.file, "\n"))
cat(paste0("Bim file is: ", bim.file, "\n"))

if (verbose) {
    cat("Continue (Y/N)? ")
    resp <- readLines(con="stdin", 1)

    if (toupper(resp)!='Y') {
        q()
    }
}

snp.info <- fread(bim.file)
setnames(snp.info, c('chr', 'snpid', 'cm', 'pos', 'allele1', 'allele2'))
setkey(snp.info, snpid)

## cat("\n\nRandomly select 20 SNPs to check.\n")

n.skip <- 0
tmp.file <- file(genotype.file)
open(tmp.file)
while (TRUE) {
    tmp.line<-readLines(tmp.file, n=1)
    if (substr(tmp.line, 1, 1)!='#') {
        break
    }
    n.skip <- n.skip+1    
}
close(tmp.file)

genotype.data <- read.delim(genotype.file, stringsAsFactors=FALSE, header=FALSE, nrow=3, skip=n.skip+1)
genotype.colnames <- scan(genotype.file, what='character', nline=1, quiet=TRUE, skip=n.skip)
if (ncol(genotype.data)==(length(genotype.colnames))) {
    names(genotype.data) <-  genotype.colnames
} else {
    if (ncol(genotype.data)==(length(genotype.colnames)+1)) {
        names(genotype.data) <- c('snpid', genotype.colnames)
    } else {
        print(genotype.colnames)
        stop("Colnames do not fit...")
    }
}

## one character one byte
title.offset <- (sum(nchar(names(genotype.data))) + ncol(genotype.data))
row.size <- sum(nchar(unlist(genotype.data[1,]))) + ncol(genotype.data)

genotype.filesize <- file.info(genotype.file)$size[1]
sample.pos <- sample(title.offset:(genotype.filesize-row.size*2), 25)

genotype.bf <- file(genotype.file, 'rb')

for (curr.pos in sample.pos) {
    seek(genotype.bf, curr.pos)

    temp <- unlist(strsplit(readChar(genotype.bf, row.size * 3), "\n"))
    if (length(temp) > 2) {
        genotype.data <- rbind(genotype.data, unlist(strsplit(temp[2], '\t')))
    }
}

close(genotype.bf)

row.names(genotype.data) <- genotype.data[,1]
genotype.data <- genotype.data[,-1]

genotype.data <- subset(genotype.data, row.names(genotype.data) %in% snp.info$snpid)

genotype.data <- genotype.data[order(row.names(genotype.data)),]

cat(paste0("Selected ", nrow(genotype.data), " SNPs for checking.\n\n"))

cat("Retrieving data from plink bed file...\n")

plink.genotype.orig <- read.plink(bed.file, select.snps=row.names(genotype.data))

plink.genotype <- t(as(plink.genotype.orig$genotype, 'character'))

plink.genotype <- plink.genotype[order(row.names(plink.genotype)),]

cat("\nCheck all snps are retrieved....\n")
all(row.names(plink.genotype)==row.names(genotype.data))

cat("\nCheck all samples are retrieved....\n")
all(colnames(plink.genotype)==names(genotype.data))

for (i in 1:nrow(genotype.data)) {
    print(c(row.names(plink.genotype)[i], row.names(genotype.data)[i]))
    print(snp.info[row.names(plink.genotype)[i],])
    print(table(unlist(plink.genotype[i,]), unlist(genotype.data[i,])))
}


q()

