plink_utilities
===============

Some scripts to help with PLINK 


*Rscript ped_maker.R genotype.txt map.file [--noprompt]*

The above R script takes a matrix format genotype file, (SNP by subject) and generate PLINK ped file. Assume MAP file has been generated. [--noprompt] will not prompt for any interactive responses.

*Rscript bed_checker.R genotype.txt bed.file [--noprompt]*

The above R script randomly select ~20 SNPs from the original genotype.txt file and cross-tabulate with generated bed file.






