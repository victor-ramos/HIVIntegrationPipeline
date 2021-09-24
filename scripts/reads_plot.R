library(ggplot2)
library(reshape2)
library(ggrepel)
library(svglite)
library(argparser)

p <- arg_parser("Integration Analysis")

p <- add_argument(p, "--input", help="stats file")

p <- add_argument(p, "--output", help="output to export the image", default="output.txt")

argv <- parse_args(p)

#tbl = read.table("~/D844K-stats.txt", header = T)
analysis.name = gsub(".*\\/(\\S+)-stats.txt","\\1", argv$input)
tbl = read.table(argv$input, header = T)
tbl.melt = melt(tbl)

plot = ggplot( data = tbl.melt, aes( x= variable, y = value, group = sample, color = sample) ) + 
  facet_wrap(. ~ sample) +
  geom_line() +
  geom_point() + 
  geom_text_repel(aes(label = value), color = "black") + 
  ggtitle(analysis.name) +
  ylab("# of reads") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank() )

#svglite(file = "~/D844K-Reads-Plot.svg", width = 12)
svglite(file = paste0(argv$output, "/", analysis.name, "-Reads-Plot.svg"), height = 10, width = 15)
print(plot)
dev.off()
