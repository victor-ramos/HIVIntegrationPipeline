#!/bin/sh

raw_data_folder=$1
output_dir=$2

for i in $(ls $raw_data_folder)
	do 
		echo $i

		mkdir -p $output_dir/hiv_integration_analysis/$i/data

		for k in $(ls $output_dir/IntegrationAnalysisPipeline)
			do
				ln -s $output_dir/IntegrationAnalysisPipeline/$k $output_dir/hiv_integration_analysis/$i 
			done

		
		for j in $(find $raw_data_folder/$i -iname "*.fastq.gz")
			do
				ln -s $j $output_dir/hiv_integration_analysis/$i/data
			done

       	done

for i in $(ls $HOME/hiv_integration_analysis/); do cd $HOME; cmd="$PWD/hiv_integration_analysis/$i/"; cd $cmd; ./scripts/run_integration_pipeline.sh; cd $HOME; done
