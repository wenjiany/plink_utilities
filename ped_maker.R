
library(data.table)

verbose <- TRUE

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

cat(paste0("Genotype file is: ", genotype.file, "\n"))
cat(paste0("Map file is: ", map.file, ".bed\n"))

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

map.data <- fread(map.file, stringsAsFactors=FALSE, header=FALSE)
setnames(map.data, c('chr', 'snpid', 'mm', 'pos'))

cat("\nReading in genotypes data. May take a long time for large dataset.\n")

genotype.data <- fread(genotype.file, stringsAsFactors=FALSE, skip=1, header=FALSE)

genotype.colnames <- scan(genotype.file, what='character', nline=1, quiet=TRUE)
if (ncol(genotype.data)==(length(genotype.colnames))) {
    setnames(genotype.data, genotype.colnames)
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

setkey(genotype.data, snpid)
genotype.data <- genotype.data[map.data$snpid,]

### possible to change to ACGT based on annotation
genotype.lookup <- c('N N', 'A A', 'A B', 'B B')

cat("\n\nGenerating target ped file: ", ped.file, "\n")

zz <- file(ped.file, "w")

for (i in 2:ncol(genotype.data)) {
    familyid <- gsub('^X', '', names(genotype.data)[i])
    subjectid <- familyid
    ## curr.sex <- inferred.sex[i]
    curr.sex <- -9
    
    curr.genotype <- unlist(genotype.data[, i, with=FALSE])
    curr.genotype <- gsub('\\/', '', curr.genotype)
    curr.genotype[!curr.genotype %in% c('AA', 'AB', 'BB')] <- 'NoCall'
                                   
    new.genotype <- genotype.lookup[as.numeric(factor(curr.genotype, levels=c('NoCall', 'AA', 'AB', 'BB')))]
    
    cat(paste(c(familyid, subjectid, 0, 0, curr.sex, unlist(new.genotype)), collapse=" "), file=zz)
    cat('\n', file=zz)
    cat('.')
}

close(zz)

cat("\n\nPed file ", ped.file, " completed\n")

map.file <- gsub('.ped$', '.map', ped.file)
cat("\nWrite map file to ", map.file, "\n")
write.table(map.data, file=map.file, sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)

cat("1) generate bed file\n")
cat(paste0("plink --noweb --no-pheno --file ", genotype.file, " --out b_", genotype.file, " --make-bed --missing-genotype N\n"))

cat("2) check consistency between bed file and original text file\n")
cat(paste0("Rscript bed_checker.R ", genotype.file, " b_", genotype.file, "\n"))

cat("3) basic SNP level summaries using bed\n")
cat(paste0("plink --noweb --bfile b_", genotype.file, " --freq \n\n"))    

q()

