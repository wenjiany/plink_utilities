
library(data.table)

if (!exists("genotype.file") || !exists("map.file")) {
    
  args.orig <- commandArgs(trailingOnly=TRUE)
  
  if ('--noprompt' %in% args.orig) {
      verbose <- FALSE
      args.orig <- setdiff(args.orig, '--noprompt')
  }
  
  if (length(args.orig) != 2) {
      stop("Usage: Rscript ped_maker.R genotype.file map.file")
  }
  
  genotype.file <- args.orig[1]
  map.file <- args.orig[2]
  
}

if (!exists("verbose")) {verbose <- TRUE}

cat(paste0("Genotype file is: ", genotype.file, "\n"))
cat(paste0("Map file is: ", map.file, "\n"))

if (verbose) {
    cat("Continue (Y/N)? ")
    resp <- readLines(con="stdin", 1)

    if (toupper(resp)!='Y') {
        q()
    }
}

cat("\n")
ped.file <- paste0(basename(genotype.file), '.ped')

if (file.exists(ped.file)) {
    if (verbose) {
        cat("ped file exists. Continue (Y/N)? ")
        resp <- readLines(con="stdin", 1)
        
        if (toupper(resp)!='Y') {
            q()
        }
    }
}

cat("Reading map file...\n")

map.data <- fread(map.file, stringsAsFactors=FALSE, header=FALSE)

if (ncol(map.data)==6) {
    setnames(map.data, c('chr', 'snpid', 'mm', 'pos', 'allele1', 'allele2'))
} else {
    setnames(map.data, c('chr', 'snpid', 'mm', 'pos'))
    map.data[, allele1:='A']
    map.data[, allele2:='B']    
}

cat("\nReading in genotypes data. May take a long time for large dataset.\n")

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

genotype.data <- fread(genotype.file, stringsAsFactors=FALSE, skip=n.skip+1, header=FALSE)

genotype.colnames <- scan(genotype.file, what='character', nline=1, quiet=TRUE, skip=n.skip)

if (ncol(genotype.data)==(length(genotype.colnames))) {
    setnames(genotype.data, genotype.colnames)
    ### for illumina data
    if (genotype.colnames[2]=='Name') {
      genotype.data <- genotype.data[, -1, with=FALSE]
      genotype.colnames <- names(genotype.data)
    }
    setnames(genotype.data, genotype.colnames[1], 'snpid')
} else {
    if (ncol(genotype.data)==(length(genotype.colnames)+1)) {
        setnames(genotype.data, c('snpid', genotype.colnames))
    } else {
        print(genotype.colnames)
        stop("Colnames do not fit...")
    }
}

print(str(genotype.data))

cat("\n\nCheck if SNPids in genotype data are in map files\n")
print(table(genotype.data$snpid %in% map.data$snpid))

if (verbose) {
    cat("\nContinue (Y/N)? ")
    resp <- readLines(con="stdin", 1)
    
    if (toupper(resp)!='Y') {
        q()
    }
}

genotype.data <- subset(genotype.data, snpid %in% map.data$snpid)

setkey(genotype.data, snpid)
genotype.data <- genotype.data[map.data$snpid,]

### possible to change to ACGT based on annotation
genotype.lookup.mat <- cbind('N N', paste(map.data$allele1, map.data$allele1),
                             paste(map.data$allele1, map.data$allele2),
                             paste(map.data$allele2, map.data$allele2))
## genotype.lookup <- c('N N', 'A A', 'A B', 'B B')

cat("\n\nGenerating target ped file: ", ped.file, "\n")

zz <- file(ped.file, "w")

for (i in 2:ncol(genotype.data)) {
    familyid <- gsub('^X', '', names(genotype.data)[i])
    subjectid <- familyid
    ## curr.sex <- inferred.sex[i]
    curr.sex <- -9
    
    curr.genotype <- unlist(genotype.data[, i, with=FALSE])

    if (any(curr.genotype %in% c(0, 1, 2))) {
        curr.genotype <- as.numeric(curr.genotype)
        curr.genotype <- curr.genotype + 2   ## AA,AB,BB==2,3,4
        curr.genotype[is.na(curr.genotype)] <- 1  ## Nocall = 1
    } else {
        curr.genotype <- gsub('\\/', '', curr.genotype)
        curr.genotype[!curr.genotype %in% c('AA', 'AB', 'BB')] <- 'NoCall'

        curr.genotype <- as.numeric(factor(curr.genotype, levels=c('NoCall', 'AA', 'AB', 'BB')))
    }
    
##   new.genotype <- genotype.lookup[as.numeric(factor(curr.genotype, levels=c('NoCall', 'AA', 'AB', 'BB')))]

    new.genotype <- genotype.lookup.mat[cbind(1:nrow(genotype.lookup.mat), curr.genotype)]
    
    cat(paste(c(familyid, subjectid, 0, 0, curr.sex, unlist(new.genotype)), collapse=" "), file=zz)
    cat('\n', file=zz)
    cat('.')
}

close(zz)

cat("\n\nPed file ", ped.file, " completed\n")

map.file <- gsub('.ped$', '.map', ped.file)
cat("\nWrite map file to ", map.file, "\n")
write.table(map.data[, 1:4, with=FALSE], file=map.file, sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)

plink_filename <- gsub('.txt$', '', basenames(genotype.file))

cat("1) generate bed file\n")
cat(paste0("plink --no-pheno --file ", plink_filename, " --out ", plink_filename, " --make-bed --missing-genotype N\n"))

cat("2) check consistency between bed file and original text file\n")
cat(paste0("Rscript bed_checker.R ", genotype.file, " ", plink_filename, "\n"))

cat("3) basic SNP level summaries using bed\n")
cat(paste0("plink --bfile ", plink_filename, " --freq \n\n"))

## q()

