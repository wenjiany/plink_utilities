plink_utilities
===============

Some scripts to help with PLINK 

* Generate ped file

Rscript ped_maker.R genotype.txt map.file [--noprompt]*

The above R script takes a matrix format genotype file, (SNP by subject) and generate PLINK ped file. Assume MAP file has been generated. [--noprompt] will not prompt for any interactive responses.

the generated ped file does not have phenotype column, therefore specify "--no-pheno" in later plink command arguments.

Note that if map file don't have allele information, the PLINK file will be generated as "A/B"; if the map file have allele information in the last two columns, the script will try to use that.

* After ped file has been generated; note the --no-pheno argument.

plink --noweb --no-pheno --file genotype.file  --out b_genotype.file --make-bed --missing-genotype N

* check generated bed file

Rscript bed_checker.R genotype.txt bed.file [--noprompt]*

The above R script randomly select ~20 SNPs from the original genotype.txt file and cross-tabulate with generated bed file.

The original genotype file is SNP by subject and coded as 'AA', 'AB' or 'BB'. ('A/A', 'A/B', or 'B/B')



