library(argparser)
library(dplyr)
library(digest)

################################ Getting argument ################################

p <- arg_parser("A text file modifying program")
p <- add_argument(p, "--bed_file", help="input file")
p <- add_argument(p, "--output", help="input file")
p <- add_argument(p, "--readsidout", help="input file")

argv = parse_args(p)

bed.file.path = argv$bed_file
output.file = argv$output
reads.id.out = argv$readsidout

############################# Filtering BED function #############################

filter.bed.file = function (bed.file.path, output.file, reads.id.out) {
  
  bed.file = read.table(file = bed.file.path, stringsAsFactors = FALSE)
  colnames(bed.file) = c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
  bed.file$score = as.numeric(bed.file$score)
  bed.file$blockSizes = as.character(bed.file$blockSizes)
  bed.file$filter = paste0(bed.file$chrom, ",", bed.file$strand, ",", bed.file$start, ",", bed.file$end)
  
  bed.file.sorted = bed.file[with(bed.file, order(chrom, -start, -end, name )), ]
  
  bed.file.sorted.filtered = bed.file.sorted %>% group_by( filter ) %>% mutate(name=paste(name,collapse = ","),score = sum(score)) %>% slice(1)
  
  bed.file.sorted.filtered.score = subset(bed.file.sorted.filtered, score >= 2)
  
  rescue.reads = data.frame(reads_id = gsub("(\\w+)?[\\,].*", "\\1", bed.file.sorted.filtered.score$name ))
  write.table(x = rescue.reads, file = reads.id.out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  bed.file.sorted.filtered.score$name = sapply(bed.file.sorted.filtered.score$name, FUN = digest, algo="md5" )
  
  bed.file.sorted.filtered.score = bed.file.sorted.filtered.score %>% 
    mutate(
      itemRgb = case_when(
      score == 2 ~ "255,0,0",
      score %in% c(3,4) ~ "240,240,0",
      score %in% c(5,6) ~ "240,140,0",
      score %in% c(7,8) ~ "0,248,0",
      score >= 9 ~ "0,140,0"
    ))
  
  print(paste0("NROW INPUT: ", nrow(bed.file.sorted)))
  
  bed.file = bed.file[!duplicated(bed.file$filter), ]
  bed.file$filter = NULL
  
  print(paste0("NROW OUTPUT: ", nrow(bed.file.sorted.filtered)))
  
  output.file.score = gsub("\\.bed", "\\_score\\.bed", output.file) 
  print(output.file.score)
  write.table(x = bed.file.sorted.filtered[, 1:12], file = output.file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(x = bed.file.sorted.filtered.score[, 1:12], file = output.file.score, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
}

#################################### Execution ##################################

filter.bed.file(bed.file.path = bed.file.path, output.file = output.file, reads.id.out = reads.id.out)
