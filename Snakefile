import os, re, glob

CWD = os.getcwd()

configfile: CWD + '/config.yml'

DATA = config['DATADIR']
OUTPUTDIR = config['OUTPUTDIR']
LTR_REF = config['LTR_REF']

FASTQ_DATA = glob.glob(CWD + '/' + DATA + '/*.fastq.gz')
FASTQ_DATA_R1 = glob.glob(CWD + '/' + DATA + '/*_R1_*.fastq.gz')
FASTQ_DATA_R2 = glob.glob(CWD + '/' + DATA + '/*_R2_*.fastq.gz')

# OUTPUT_BAIT = [ CWD + '/' + OUTPUTDIR + '/pre_selected_reads/' + re.sub("_L.*", "", os.path.basename(file)) + '/' + re.sub(r'(\S+)(\_L\d+\_.*)',r'\1_bait\2' , os.path.basename(file))  for file in FASTQ_DATA_R1 ]
# OUTPUT_TARGET = [ CWD + '/' + OUTPUTDIR + '/pre_selected_reads/' + re.sub("_L.*", "", os.path.basename(file)) + '/' + re.sub(r'(\S+)(\_L\d+\_.*)',r'\1_target\2' , os.path.basename(file))  for file in FASTQ_DATA_R2 ]


FINAL_OUTPUT = [ CWD + '/' + OUTPUTDIR + '/fragment_count/' + re.sub("_L.*", "", os.path.basename(file)) + '/' + re.sub(r'_R[12].*?(\d+)\.fastq\.gz',r'_\1_frag_count.txt' , os.path.basename(file))  for file in FASTQ_DATA ]

RESCUE_FASTQ_R1_OUTPUT = [ CWD + '/' + OUTPUTDIR + '/rescued_reads/' + re.sub("_L.*", "", os.path.basename(file)) + '/' + re.sub(r'(\w+)_L(\d+)_R1_(\w+)\.fastq\.gz',r'\1_bait_rescue_L\2_R1_\3.fq' , os.path.basename(file))  for file in FASTQ_DATA_R1 ]
RESCUE_FASTQ_R2_OUTPUT = [ CWD + '/' + OUTPUTDIR + '/rescued_reads/' + re.sub("_L.*", "", os.path.basename(file)) + '/' + re.sub(r'(\w+)_L(\d+)_R2_(\w+)\.fastq\.gz',r'\1_target_rescue_L\2_R2_\3.fq' , os.path.basename(file))  for file in FASTQ_DATA_R2 ]


rule run_all:
    input: [FINAL_OUTPUT, RESCUE_FASTQ_R1_OUTPUT, RESCUE_FASTQ_R2_OUTPUT]
    # input: [OUTPUT_BAIT, OUTPUT_TARGET]


rule select_reads_cutadapt:
    input:
        expand(CWD + '/' + DATA + '/{{sample}}_L{{lane}}_R{pair}_{{suffix}}.fastq.gz', pair=['1','2'])
    output:
        bait_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_bait_L{lane}_R1_{suffix}.fq",
        target_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_target_L{lane}_R2_{suffix}.fq",
        bait_fastq_aux = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_bait_L{lane}_R1_{suffix}_aux.fastq.gz",
        target_fastq_aux = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_target_L{lane}_R2_{suffix}_aux.fastq.gz",
        bait_fastq_sorted = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_bait_L{lane}_R1_{suffix}_sorted.fastq",
        target_fastq_sorted = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_target_L{lane}_R2_{suffix}_sorted.fastq",
        unkown_reads_r1_ref = expand(CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{{sample}}/{{sample}}_unknown_R1_ref_L{{lane}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2']),
        unkown_reads_r2_ref = expand(CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{{sample}}/{{sample}}_unknown_R2_ref_L{{lane}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2'])
    log:
        bait_r1_r1_ref = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_L{lane}_{suffix}_bait_R1_R1_ref.log",
        bait_unknown_r2_r1_ref = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_L{lane}_{suffix}_bait_unknown_R2_R1_ref.log",
        linker_r2_r2_ref = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_L{lane}_{suffix}_linker_R2_R2_ref.log",
        linker_unknown_r1_r2_ref = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_L{lane}_{suffix}_linker_unknown_R1_R2_ref.log"
    params:
        dir = CWD+"/"+OUTPUTDIR+"/pre_selected_reads",
        bait_target_common_reads_list = CWD + "/" + OUTPUTDIR + "/pre_selected_reads/{sample}/bait_target_commom_reads.list"
    shell:
        """
            FIRST_CUTADAPT_OUTPUT_R1='{params.dir}/{wildcards.sample}/{wildcards.sample}_bait_R1_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz'
            FIRST_CUTADAPT_OUTPUT_R2='{params.dir}/{wildcards.sample}/{wildcards.sample}_bait_R1_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz'

            #Subseting reads that start with the LTR in R1
            cutadapt --action=lowercase --no-indels -e 0 -g "bait=XAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA;min_overlap=34" \
               -o {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_R1_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz \
               -p {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_R1_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz \
               {input[0]} \
               {input[1]} &> {log.bait_r1_r1_ref}


            SECOND_CUTADAPT_OUTPUT_R1='{params.dir}/{wildcards.sample}/{wildcards.sample}_bait_unknown_R2_R1_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz'
            SECOND_CUTADAPT_OUTPUT_R2='{params.dir}/{wildcards.sample}/{wildcards.sample}_bait_unknown_R2_R1_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz'

            #Subseting reads that start with the LTR in unkown R2 (R1 ref)
            cutadapt --action=lowercase --no-indels -e 0 -g "bait=XAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA;min_overlap=34" \
               -o {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_unknown_R2_R1_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz \
               -p {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_unknown_R2_R1_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz \
               {output.unkown_reads_r1_ref[1]} \
               {output.unkown_reads_r1_ref[0]} &> {log.bait_unknown_r2_r1_ref}



            THIRD_CUTADAPT_OUTPUT_R1='{params.dir}/{wildcards.sample}/{wildcards.sample}_linker_R2_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz'
            THIRD_CUTADAPT_OUTPUT_R2='{params.dir}/{wildcards.sample}/{wildcards.sample}_linker_R2_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz'

            #Subseting reads that start with the Linker in R2
            cutadapt --action=lowercase --no-indels -e 0 -g "linker=XGCAGCGGATAACAATTTCACACAGGACGTACTGTGGCGCGCCT;min_overlap=43" \
                -o {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_R2_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz \
                -p {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_R2_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz \
                {input[1]} \
                {input[0]} &> {log.linker_r2_r2_ref}


            FOURTH_CUTADAPT_OUTPUT_R1='{params.dir}/{wildcards.sample}/{wildcards.sample}_linker_unknown_R1_R2_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz'
            FOURTH_CUTADAPT_OUTPUT_R2='{params.dir}/{wildcards.sample}/{wildcards.sample}_linker_unknown_R1_R2_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz'

            #Subseting reads that start with the Linker in unkown R1 (R2 ref)
            cutadapt --action=lowercase --no-indels -e 0 -g "linker=XGCAGCGGATAACAATTTCACACAGGACGTACTGTGGCGCGCCT;min_overlap=43" \
               -o {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_unknown_R1_R2_ref_L{wildcards.lane}_R1_{wildcards.suffix}.fastq.gz \
               -p {params.dir}/{wildcards.sample}/{wildcards.sample}_{{name}}_unknown_R1_R2_ref_L{wildcards.lane}_R2_{wildcards.suffix}.fastq.gz \
               {output.unkown_reads_r2_ref[1]} \
               {output.unkown_reads_r2_ref[0]} &> {log.linker_unknown_r1_r2_ref}


            if [[ -e $FIRST_CUTADAPT_OUTPUT_R1 && -e $SECOND_CUTADAPT_OUTPUT_R1 ]];then
                cat <( gunzip -c $FIRST_CUTADAPT_OUTPUT_R1 ) \
                    <( gunzip -c $SECOND_CUTADAPT_OUTPUT_R1 ) \
                    | gzip - > {output.bait_fastq_aux}
            else
                cp $FIRST_CUTADAPT_OUTPUT_R1 {output.bait_fastq_aux}
            fi

            if [[ -e $THIRD_CUTADAPT_OUTPUT_R1 && -e $FOURTH_CUTADAPT_OUTPUT_R1 ]];then
                cat <( gunzip -c $THIRD_CUTADAPT_OUTPUT_R1 ) \
                    <( gunzip -c $FOURTH_CUTADAPT_OUTPUT_R1 ) \
                    | gzip - > {output.target_fastq_aux}
            else
                cp $THIRD_CUTADAPT_OUTPUT_R1 {output.target_fastq_aux}
            fi


            #cat <(gunzip -c {output.bait_fastq_aux} | grep '@M') <(gunzip -c {output.target_fastq_aux} ) | perl -lane 'print $1 if $_ =~ m/(\S+)\s(\S+)/ ' | sort | uniq -c | perl -lane 'print $F[1] if $F[0] > 1 ' > {params.bait_target_common_reads_list}
            cat <( bp_seqconvert --from fastq --to tab < <( gunzip -c {output.bait_fastq_aux} ) | cut -f 1  ) <( bp_seqconvert --from fastq --to tab < <( gunzip -c {output.target_fastq_aux} ) | cut -f 1 ) | sort | uniq -c | perl -lane 'print $F[1] if $F[0] > 1 ' > {params.bait_target_common_reads_list}
            #sed -i 's/\@M/M/g' {params.bait_target_common_reads_list}

            seqtk subseq {output.bait_fastq_aux} {params.bait_target_common_reads_list} > {output.bait_fastq_sorted}
            seqtk subseq {output.target_fastq_aux} {params.bait_target_common_reads_list} > {output.target_fastq_sorted}

            fastq-sort --id {output.bait_fastq_sorted} > {output.bait_fastq}
            fastq-sort --id {output.target_fastq_sorted} > {output.target_fastq}

        """


rule ltr_alignment:
    input:
        bait_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_bait_L{lane}_R1_{suffix}.fq",
        target_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_target_L{lane}_R2_{suffix}.fq"
    output: 
        sam = CWD + '/' + OUTPUTDIR + '/smalt_results/{sample}/{sample}_L{lane}_bait_target_{suffix}.sam',
    threads: 3
    resources:
        mem_gb = 3
    params:
        ref = LTR_REF,
        match = config['MATCH'],
        subst=config['SUBST'],
        gapopen=config['GAPOPEN'],
        gapext=config['GAPEXT']
    log:
        smalt_log = CWD + '/' + OUTPUTDIR + '/smalt_results/{sample}/{sample}_L{lane}_R12_{suffix}_smalt.log',
    shell:
        """
           smalt map -f sam:x -S match={params.match},subst={params.subst},gapopen={params.gapopen},gapext={params.gapext} {params.ref} {input.bait_fastq} {input.target_fastq} 1> {output} 2> {log}
        """


rule get_real_baits:
    input: 
        sam = CWD + '/' + OUTPUTDIR + '/smalt_results/{sample}/{sample}_L{lane}_bait_target_{suffix}.sam'
    output:
        ltr_reads_to_recover = CWD + '/' + OUTPUTDIR + '/recover_ltr_reads/{sample}/{sample}_L{lane}_bait_target_{suffix}_bait.list'
    shell:
        """
            samtools view -q 20 -F 16 {input} | perl -lane 'print $F[0] if $F[5] =~ /^34M/' > {output}
        """


rule get_bait_and_target_subset:
    input:
        #expand(CWD + '/' + OUTPUTDIR + '/recover_ltr_reads/{{sample}}/{{sample}}_L{{lane}}_R{pair}_{{suffix}}_bait.list', pair=['1','2'])
        ltr_reads_to_recover = CWD + '/' + OUTPUTDIR + '/recover_ltr_reads/{sample}/{sample}_L{lane}_bait_target_{suffix}_bait.list'
    output:
        expand(CWD+"/"+OUTPUTDIR+"/bait_target_fastq/{{sample}}/{{sample}}_L{{lane}}_{pair}_{{suffix}}.fq",pair=['bait','target'])
    params:
        bait_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_bait_L{lane}_R1_{suffix}.fq",
        target_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_target_L{lane}_R2_{suffix}.fq"


    shell:
        """
            seqtk subseq {params.bait_fastq} {input.ltr_reads_to_recover} > {output[0]}
            seqtk subseq {params.target_fastq} {input.ltr_reads_to_recover} > {output[1]}

        """


rule append_bait_and_target:
    input:
        bait = CWD+"/"+OUTPUTDIR+"/bait_target_fastq/{sample}/{sample}_L{lane}_bait_{suffix}.fq",
        target = CWD+"/"+OUTPUTDIR+"/bait_target_fastq/{sample}/{sample}_L{lane}_target_{suffix}.fq"
    output:
        expand(CWD+"/"+OUTPUTDIR+"/append_bait_target_fastq/{{sample}}/{{sample}}_L{{lane}}_{pair}_{{suffix}}.fq",pair=['bait','target'])
    shell:
        """
            cat {input.bait} > {output[0]}
            cat {input.target} > {output[1]}
            perl -i -lane 'if ($_ =~ /^(@\S+)?(\s[12].*)/){{print $1."_bait".$2; }} else {{ print }} ' {output[0]}
            perl -i -lane 'if ($_ =~ /^(@\S+)?(\s[12].*)/){{print $1."_target".$2; }} else {{ print }} ' {output[1]}
        """


rule ltr_and_linker_trimming:
    input:
        expand(CWD+"/"+OUTPUTDIR+"/append_bait_target_fastq/{{sample}}/{{sample}}_L{{lane}}_{pair}_{{suffix}}.fq",pair=['bait','target'])
    output:
        bait_fastq = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_bait_{suffix}_aux.fq",
        target_fastq = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_target_{suffix}_aux.fq",
        bait_fastq_aux = temp(CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_bait_{suffix}_val_1_aux.fq"),
        target_fastq_aux = temp(CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_target_{suffix}_val_2_aux.fq"),
        bait_fastq_final = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_bait_{suffix}_val_1.fq",
        target_fastq_final = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_target_{suffix}_val_2.fq"
    params:
        ltr_34_fasta = config['LTR_FASTA'],
        ltr_34_rev_comp_fasta = config['LTR_REV_COMP_FASTA'],
        linker_fasta = config['LINKER_FASTA'],
        linker_rev_comp_fasta = config['LINKER_REV_COMP_FASTA']
    log:
        ltr_34 = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_{suffix}_ltr_34_trim.log",
        linker_rev_comp = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_{suffix}_linker_rev_comp_trim.log",
        linker = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_{suffix}_linker_trim.log",
        ltr_34_rev_comp = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_{suffix}_ltr_34_rev_comp_trim.log"
    threads: 2
    resources:
        mem_mb=10000
    shell:
        """

            bbduk.sh in={input[0]} ref={params.ltr_34_fasta} k=15 ktrim=l minlength=0 out={output.bait_fastq} &> {log.ltr_34}
            bbduk.sh in={output.bait_fastq} ref={params.linker_rev_comp_fasta} k=15 ktrim=r qtrim=l trimq=20 minlength=0 out={output.bait_fastq_aux} &> {log.linker_rev_comp}
            bbduk.sh in={input[1]} ref={params.linker_fasta} k=15 ktrim=l minlength=0 out={output.target_fastq} &> {log.linker}
            bbduk.sh in={output.target_fastq} ref={params.ltr_34_rev_comp_fasta} k=15 ktrim=r qtrim=l trimq=20 minlength=0 out={output.target_fastq_aux} &> {log.ltr_34_rev_comp}

            fastq-sort --id {output.bait_fastq_aux} > {output.bait_fastq_final}
            fastq-sort --id {output.target_fastq_aux} > {output.target_fastq_final} 

        """


rule align_to_human_and_hiv_genome:
    input:
        input_bait = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_bait_{suffix}_val_1.fq",
        input_target = CWD+"/"+OUTPUTDIR+"/bbduk_trim/{sample}/{sample}_L{lane}_target_{suffix}_val_2.fq"
    output:
        #bed=CWD+"/"+OUTPUTDIR+"/reads_alignment_bed/{sample}/{sample}_L{lane}_{suffix}_bait_target.bed"
        bam=CWD+"/"+OUTPUTDIR+"/reads_alignment/{sample}/{sample}_L{lane}_{suffix}_bait_target.bam",
        bam_hiv=CWD+"/"+OUTPUTDIR+"/reads_alignment/{sample}/{sample}_L{lane}_{suffix}_bait_target_hiv.bam",
        hiv_bait_fastq = CWD+"/"+OUTPUTDIR+"/reads_alignment/{sample}/{sample}_L{lane}_bait_{suffix}_hiv.fq",
        hiv_target_fastq = CWD+"/"+OUTPUTDIR+"/reads_alignment/{sample}/{sample}_L{lane}_target_{suffix}_hiv.fq"
    params: 
        genome_index = config['HUMAN_REF'],
        hiv_genome_index = config['HXB2_REF'],
        unmapped_bait = CWD + "/" + OUTPUTDIR + "/reads_alignment/{sample}/{sample}_L{lane}_{suffix}_unmapped_bait.list",
        unmapped_target = CWD + "/" + OUTPUTDIR + "/reads_alignment/{sample}/{sample}_L{lane}_{suffix}_unmapped_target.list"
    threads: 24
    shell:
        """
            smalt map -x -r -1 -f sam:x -n {threads} {params.genome_index} {input.input_bait} {input.input_target} | samtools view -Sb > {output.bam}

            samtools view -F 2 {output.bam} | cut -f1 | grep 'bait' > {params.unmapped_bait}
            samtools view -F 2 {output.bam} | cut -f1 | grep 'target' > {params.unmapped_target}

            seqtk subseq {input.input_bait} {params.unmapped_bait} | fastq-sort --id > {output.hiv_bait_fastq}
            seqtk subseq {input.input_target} {params.unmapped_target} | fastq-sort --id > {output.hiv_target_fastq}

            smalt map -f sam:x  -n {threads} {params.hiv_genome_index} {output.hiv_bait_fastq} {output.hiv_target_fastq} | samtools view -Sb > {output.bam_hiv}

        """

rule filter_human_reads:
    input:
        bam=CWD+"/"+OUTPUTDIR+"/reads_alignment/{sample}/{sample}_L{lane}_{suffix}_bait_target.bam"
    output:
        header_tmp_mapped_paired_reads = temp(CWD+"/"+OUTPUTDIR+"/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_header_mapped_bait_target.sam"),
        mapped_paired_reads = CWD+"/"+OUTPUTDIR+"/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_proper_mapped_bait_target.sam",
        tmp_filtered_mapped_paired_reads = temp(CWD+"/"+OUTPUTDIR+"/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_tmp_filtered_mapped_bait_target.sam"),
        filtered_mapped_paired_reads = CWD+"/"+OUTPUTDIR+"/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_filtered_mapped_bait_target.sam",
        bed6_human_alignment = CWD + "/" + OUTPUTDIR + "/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_alignmed_reads_raw.bed",
        raw_mapped_paired_reads_tmp = temp(CWD+"/"+OUTPUTDIR+"/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_raw_mapped_bait_target_tmp.sam")
    params:
        insertmax = 2000
    shell:
        """
            samtools view -H  {input.bam} > {output.header_tmp_mapped_paired_reads}
            samtools view {input.bam} > {output.raw_mapped_paired_reads_tmp}

            bamToBed -i {input.bam} > {output.bed6_human_alignment}

            perl scripts/getProperPairsAndOverlappingReads.pl get_proper_pairs_and_overlapping_reads \
            --insertmax {params.insertmax} \
            --bed6 {output.bed6_human_alignment} \
            --sam_file {output.raw_mapped_paired_reads_tmp} > {output.mapped_paired_reads}

            perl scripts/FilterHumanAlignment.pl filter_human_alignment --input_file {output.mapped_paired_reads} > {output.tmp_filtered_mapped_paired_reads}

            cat {output.header_tmp_mapped_paired_reads} {output.tmp_filtered_mapped_paired_reads} > {output.filtered_mapped_paired_reads}
        """

rule get_bed_file_and_hotspots:
    input:
        filtered_mapped_filtered_reads = CWD+"/"+OUTPUTDIR+"/filtered_human_reads_alignment/{sample}/{sample}_L{lane}_{suffix}_filtered_mapped_bait_target.sam"
    output:
        temp_bed = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_bait_target.bed",
        bed = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_bait_target_filtered_per_fragment.bed",
        bed_score = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_bait_target_filtered_per_fragment_score.bed",
        hot_spots = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_hot_spots.tsv",
        rescue_representative_reads = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_rescue_representative_reads.list"
    params:
        bed_bait_target = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_filtered_mapped_bait_target.bed"
    log:
        CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_filtered_mapped_bait_target.log"
    shell:
        """
            samtools view -u  {input.filtered_mapped_filtered_reads} | samtools sort -n  | bamToBed -i > {params.bed_bait_target}

            if [[ -s {params.bed_bait_target} ]];then
                perl scripts/bed_converter.pl bed6_to_bed12 --input_file {params.bed_bait_target} > {output.temp_bed}

                Rscript scripts/BedFilePerFragment.R --bed_file {output.temp_bed} --output {output.bed} --readsidout {output.rescue_representative_reads} > {log}

                echo "track name='{wildcards.sample}' description='{wildcards.sample}'" itemRgb="On" > {output.bed_score}.aux
                cat {output.bed_score} >> {output.bed_score}.aux

                perl -Iscripts/lib scripts/bin/BioTools.pl cluster -t --human --minority 0 --input_file {output.bed_score} > {output.hot_spots}
            else
                touch {output.temp_bed}
                touch {output.bed}
                touch {output.bed_score}
            fi

        """


rule rescue_representative_reads_from_pre_selected:
    input:
        rescue_representative_reads = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_rescue_representative_reads.list"
    output:
        rescue_bait_fastq = CWD+"/"+OUTPUTDIR+"/rescued_reads/{sample}/{sample}_bait_rescue_L{lane}_R1_{suffix}.fq",
        rescue_target_fastq = CWD+"/"+OUTPUTDIR+"/rescued_reads/{sample}/{sample}_target_rescue_L{lane}_R2_{suffix}.fq"
    params:
        bait_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_bait_L{lane}_R1_{suffix}.fq",
        target_fastq = CWD+"/"+OUTPUTDIR+"/pre_selected_reads/{sample}/{sample}_target_L{lane}_R2_{suffix}.fq"
    shell:
        """
            seqtk subseq {params.bait_fastq} {input.rescue_representative_reads} | fastq-sort --id > {output.rescue_bait_fastq}
            seqtk subseq {params.target_fastq} {input.rescue_representative_reads} | fastq-sort --id > {output.rescue_target_fastq}
        """

rule chr_insertion_count:
    input:
        bed = CWD+"/"+OUTPUTDIR+"/human_bed_and_hostspot/{sample}/{sample}_L{lane}_{suffix}_bait_target_filtered_per_fragment_score.bed"
    output:
        txt=CWD+"/"+OUTPUTDIR+"/fragment_count/{sample}/{sample}_L{lane}_{suffix}_frag_count.txt"
    shell:
        """
            if [[ -s {input.bed}  ]];then
                cut -f1 {input.bed} | sort -k1V,1 | uniq -c > {output.txt}
            else
                touch {output.txt}
            fi
        """

