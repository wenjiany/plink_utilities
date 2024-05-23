
library(data.table)

annotation_file <- 'data/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt'

ann_data <- fread(annotation_file, sep="\t")

dim(ann_data)

tail(ann_data)

ann_data$cm <- 0

map_data <- ann_data[, c("Chr", "Name", "cm", "MapInfo")]

fwrite(map_data, file="data/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.map", col.names=FALSE, sep="\t")

