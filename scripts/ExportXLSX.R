################################## Libraries ##################################
library("xlsx")
library("argparser")


################################## Functions ##################################

export.xlsx.file = function ( file, wb ) {
  
  hot.spots.aux.file = read.table( file )
  
  sample = gsub("\\S+\\/human_bed_and_hostspot\\/(\\S+)\\/\\S+", "\\1", file)
  
  hot.spots.file = hot.spots.aux.file[, c(1, 2, 3, 4, 5, 6, 9)]
  colnames(hot.spots.file) = c("chr", "start", "end", "hotspot_id", "hotspot_size", "number_of_fragments", "p-value")
  
  hot.spots.file.final <- hot.spots.file[ rep(1:nrow( hot.spots.file ), each = 2), ]
  hot.spots.file.final[ 1:nrow( hot.spots.file.final ) %% 2 == 0, ] <- ""
  
  SUB_TITLE_STYLE = CellStyle(wb) + 
    Font(wb, isBold = TRUE) + 
    Alignment(horizontal="ALIGN_CENTER")
  
  sheet <- createSheet(wb, sheetName = sample)
  
  addDataFrame(hot.spots.file.final, sheet, row.names = FALSE, colnamesStyle = SUB_TITLE_STYLE)
  
}

################################## Execution ##################################

p <- arg_parser("Export xlsx file")

p <- add_argument(p, "--input", help="hotspots file")

p <- add_argument(p, "--output", help="xlsx file output",)

argv <- parse_args(p)

hot.spots.files = list.files( argv$input, recursive = T, full.names = T, pattern = "*hot_spots.tsv" )

wb = createWorkbook(type="xlsx")
for (file in hot.spots.files) {
    if (file.size(file) > 0){
        export.xlsx.file( file = file, wb = wb )
    }
}
saveWorkbook(wb, argv$output)

