

genotype.dir <- "Z:/ResearchHome/Groups/yanggrp/projects/YanggrpCommon/common/CommonWorkflow/GDA/"

genotype.files <- list.files(genotype.dir, pattern="202.*_GDA_.*.genotypes.txt$", full.names = TRUE)

length(genotype.files)

map.file <- "data/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.map"

verbose <- FALSE

for (genotype.file in genotype.files) {
  source("ped_maker.R")
}

