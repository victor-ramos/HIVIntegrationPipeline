#!/bin/sh

export LD_LIBRARY_PATH=/ru-auth/local/home/vramos/miniconda3/envs/integration_env/lib/server

analysisname=`basename $PWD`

echo $analysisname

snakemake -j $(nproc --all) -k -s $PWD/Snakefile

Rscript $PWD/scripts/ExportXLSX.R --input $PWD/results/human_bed_and_hostspot --output $PWD/${analysisname}_hotspots.xlsx

mkdir -p $PWD/results/plots

perl $PWD/scripts/getReads.pl get_reads --fastq_dir $PWD/data --analysis_dir $PWD/results > $PWD/results/plots/$analysisname-stats.txt

Rscript $PWD/scripts/reads_plot.R --input $PWD/results/plots/$analysisname-stats.txt --output $PWD/results/plots

Rscript $PWD/scripts/plot_counts_by_chr.R --input $PWD/results/fragment_count --analysisname $analysisname --output $PWD/results/plots/

for i in $(find $PWD/results/human_bed_and_hostspot/ -iname "*aux"); do cat $i >> $PWD/integrations_count_all_genome_broswer_$analysisname.txt; done
