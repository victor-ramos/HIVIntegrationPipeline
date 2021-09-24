library("ggplot2")
library("ggpubr")
library("gtools")
library("svglite")
library(argparser)

p <- arg_parser("Integration Analysis")

p <- add_argument(p, "--input", help="input dir")
p <- add_argument(p, "--analysisname", help="analysis name")
p <- add_argument(p, "--output", help="output to export the image", default="output.txt")

argv <- parse_args(p)

get.integration.plot = function ( conc ) {
  
  integration.plot = ggplot(conc, aes(x = sample, y = count, fill = factor(chr), group=factor(sample.name) )) +
    geom_bar(position = position_dodge2(preserve = "single"), stat = "identity", color="black" ) +
    geom_text(aes(label = count), position = position_dodge2(0.9), vjust = -1, color = "black", fontface="bold") +
    coord_cartesian(expand = FALSE, ylim = c(0,350)) +
    ggtitle("Fragments count per chromosome") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          plot.title =  element_text(hjust = 0.5))
  
  return(integration.plot)
}


################################################################################

# counts.file = list.files("~/data04-scratch/analysis/Amy_Huang/D7Y8T/results/fragment_count", recursive = TRUE, full.names = TRUE)
counts.file = list.files(argv$input, recursive = TRUE, full.names = TRUE)
counts.file = grep("txt", counts.file, value = TRUE)
conc = NULL

for (file in counts.file) {
  
  sample.name = gsub("(\\S+)_L(\\S+)_(\\S+)_frag_count.txt", "\\1", basename(file), perl = T)
  
  if ( file.size( file ) > 0 ) {
    aux = read.table(file )
    colnames(aux) = c("count","chr")
    aux$sample = sample.name
    aux$chr = factor(aux$chr, levels = sort(aux$chr) )
    
    conc = rbind(conc, aux)  
  }
}

conc$chr = factor(conc$chr, levels = unique(mixedsort(as.character(conc$chr))))

integration.plot = get.integration.plot(conc)


#svglite(file = "~/integration_plot_D7Y8T.svg", width = 6, height = 7)
svglite(file = paste0(argv$output, "integration_plot_",argv$analysisname,".svg"), width = 8, height = 10)
print(integration.plot)
dev.off()
