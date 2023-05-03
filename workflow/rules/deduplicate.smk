rule move_umi_to_rx_tag:
        input:
                bam_sorted = config['resultsdir'] + "/results/4_alignment/bam_filt/{sample}_filt.bam"
        output:
                bam_sorted_UMIrx = temp(config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.bam")
        log:
                bamtag_log = config['resultsdir'] + "/logs/5_dedup/{sample}.rxbamtag.log"
        conda:
                "envs/twist_target.yaml"
        shell:
                """
                srslyumi-bamtag \
                        --binary \
                        -o {output.bam_sorted_UMIrx} \
                        --take-fragment 0 \
                        {input.bam_sorted} \
                2> {log.bamtag_log}
                """

rule index_rxbam:
        input:
                rx_bam = config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.bam"
        output:
                config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.bam.bai"
        log:
                config['resultsdir'] + "/logs/5_dedup/{sample}.index_rxbam.log"
        threads: index_threads
        shell:
                """
                samtools index -@ {threads} {input.rx_bam}
                """

rule group_umi_umitools:
        input:
                bam_filt_UMIrx = config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.bam",
                bam_filt_UMIrx_index = config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.bam.bai"
        output:
                bam_filt_UMIgrouped = config['resultsdir'] + "/results/5_dedup/{sample}_filt.umi.bam",
                bam_filt_UMItsv = config['resultsdir'] + "/results/5_dedup/{sample}_filt.umitools.tsv",
                bam_filt_UMIgrouped_sort = config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort.umi.bam",
                bam_filt_UMIgrouped_sort_coord = config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort_coord.umi.bam"
        log:
                log_umitools = config['resultsdir'] + "/logs/5_dedup/{sample}.umitools_group.log",
                error = config['resultsdir'] + "/logs/5_dedup/{sample}.umitools_group.err",
                log_sort = config['resultsdir'] + "/logs/5_dedup/{sample}.sort.log",
                log_sort_coord = config['resultsdir'] + "/logs/5_dedup/{sample}.sort.log"
        params:
                min_map_q=0,
                random_seed=42
        conda:
                "envs/twist_target.yaml"
        threads: dedup_threads
        shell:
                """
                umi_tools group \
                        --output-bam \
                        --stdin={input.bam_filt_UMIrx} \
                        --stdout={output.bam_filt_UMIgrouped} \
                        --group-out={output.bam_filt_UMItsv} \
                        --log={log.log_umitools} \
                        --error={log.error} \
                        --extract-umi-method=tag \
                        --umi-tag=RX \
                        --method=directional \
                        --paired \
                        --unmapped-reads=use \
                        --mapping-quality={params.min_map_q} \
                        --random-seed={params.random_seed}

                samtools sort \
                        --template-coordinate \
                        -@ {threads} \
                        -o {output.bam_filt_UMIgrouped_sort} \
                        {output.bam_filt_UMIgrouped} \
                2> {log.log_sort}

                samtools sort \
                        -@ {threads} \
                        -o {output.bam_filt_UMIgrouped_sort_coord} \
                        {output.bam_filt_UMIgrouped} \
                2> {log.log_sort_coord}
                """

rule index_umi_bam:
        input:
                bam_filt_UMIgrouped_sort_coord = config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort_coord.umi.bam"
        output:
                config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort_coord.umi.bam.bai"
        log:
                config['resultsdir'] + "/logs/5_dedup/{sample}.index_umibam.log"
        threads: index_threads
        shell:
                """
                samtools index -@ {threads} {input.bam_filt_UMIgrouped_sort_coord}
                """

rule UMI_dedup:
        input:
                bam_filt_UMIgrouped_sort_coord = config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort_coord.umi.bam",
                bam_filt_UMIgrouped_sort_coord_index = config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort_coord.umi.bam.bai"
        output:
                config['resultsdir'] + "/results/5_dedup/{sample}_edit_distance.tsv",
                config['resultsdir'] + "/results/5_dedup/{sample}_per_umi_per_position.tsv",
                config['resultsdir'] + "/results/5_dedup/{sample}_per_umi.tsv",
                bam_dedup = config['resultsdir'] + "/results/5_dedup/{sample}_filt_dedup_umitools.bam"
        log:
                umitoolsdedup_log = config['resultsdir'] + "/logs/5_dedup/{sample}.dedup_umitools.log"
        params:
                grouping_method = config['umi_tools_grouping_method'],
                output_prefix = config['resultsdir'] + "/results/5_dedup/{sample}"
        conda:
                "envs/twist_target.yaml"
        shell:
                """
                umi_tools dedup \
                        --stdin={input.bam_filt_UMIgrouped_sort_coord} \
                        --output-stats={params.output_prefix} \
                        --log={log.umitoolsdedup_log} \
                        --extract-umi-method=tag \
                        --umi-tag=RX \
                        --method {params.grouping_method} \
                        --paired \
                        --chimeric-pairs=discard \
                        --unpaired-reads=discard \
                > {output.bam_dedup}
                """

rule setmateinfo:
        input:
                bam_filt_UMIrx = config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.bam"
        output:
                bam_filt_UMIrx_setmate = config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.setmate.bam"
        log:
                log_sort_query = config['resultsdir'] + "/logs/5_dedup/{sample}.sort_queryname.log",
                log_setmateinfo = config['resultsdir'] + "/logs/5_dedup/{sample}.SetMateInformation_fgbio.log"
        params:
                tmp_dir=config['tmp-dir']
        conda:
                "envs/twist_target.yaml"
        threads: sort_filt_threads
        shell:
                """
                samtools sort \
                        -n \
                        -T {params.tmp_dir} \
                        -@ {threads} \
                        -o /dev/stdout \
                        {input.bam_filt_UMIrx} \
                2> {log.log_sort_query} | \
                fgbio SetMateInformation \
                        -i /dev/stdin \
                        -o {output.bam_filt_UMIrx_setmate} \
                2> {log.log_setmateinfo}
                """

rule group_UMI_fgbio:
        input:
                bam_filt_UMIrx_setmate = config['resultsdir'] + "/results/5_dedup/{sample}_filt.rx.setmate.bam"
        output:
                bam_filt_UMIgrouped = config['resultsdir'] + "/results/5_dedup/{sample}_filt.fgbio.bam",
                histo = config['resultsdir'] + "/results/5_dedup/{sample}.fgbio.histo.tsv"
        log:
                fgbio_log = config['fastqdir'] +"/logs/5_dedup/{sample}_filt.group_fgbio.log"
        conda:
                "envs/twist_target.yaml"
        threads: dedup_threads
        shell:
                """
                fgbio GroupReadsByUmi \
                        -i {input.bam_filt_UMIrx_setmate} \
                        -o {output.bam_filt_UMIgrouped} \
                        --s adjacency \
                2> {log.fgbio_log}
                """

rule fgbio_consensus_read:
        input:
                bam_filt_UMIgrouped_sort = config['resultsdir'] + "/results/5_dedup/{sample}_filt_sort.umi.bam"
        output:
                bam_filt_cons = config['resultsdir'] + "/results/5_dedup/{sample}_filt.cons.bam"    
        log:
                fgbio_log = config['resultsdir'] +"/logs/5_dedup/{sample}.consensus.log"
        params:
                tmp_dir=config['tmp-dir'],
                reads_per_umi=config['reads_per_umi']
        resources:
                mem_mb=5000
        conda:
                "envs/twist_target.yaml"
        threads: dedup_threads
        shell:
                """
                fgbio -Xmx{resources.mem_mb}M --tmp-dir={params.tmp_dir} \
                        CallMolecularConsensusReads \
                        -i {input.bam_filt_UMIgrouped_sort} \
                        -o {output.bam_filt_cons} \
                        --tag=RX \
                        --min-reads={params.reads_per_umi} \
                        --threads={threads} \
                2> {log.fgbio_log}
                """

rule consensus_bam_to_fastq:
        input:
                bam_filt_cons = config['resultsdir'] + "/results/5_dedup/{sample}_filt.cons.bam"
        output:
                R1_cons = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.R1.fastq.gz",
                R2_cons = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.R2.fastq.gz"
        log:
        params:
        conda:
                "envs/twist_target.yaml"
        threads:
                align_threads
        shell:
                """
                samtools fastq -@ {threads} -1 {output.R1_cons} -2 {output.R2_cons} {input.bam_filt_cons}
                """

rule align_consensus_reads:
        input:
                R1_cons = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.R1.fastq.gz",
                R2_cons = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.R2.fastq.gz"      
        output:
                bam_cons = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.bam"
        log:
                align_log = config['resultsdir'] + "/logs/8_dd-cfDNA/{sample}.cons_bwameth.log",
                samb_log = config['resultsdir'] + "/logs/8_dd-cfDNA/{sample}.cons_toBam_sambamba.log" # change to toBam!!!!
        params:
                reference = config['reference_wo_ctrl']
        conda:
                "envs/twist_target.yaml"
        threads: align_threads
        shell:
                """
                bwameth.py \
                        --reference {params.reference} \
                        -t {threads} \
                        {input.R1_cons} {input.R2_cons} \
                        2> {log.align_log} | \
                sambamba view \
                        -h \
                        -t {threads} \
                        --sam-input \
                        --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 30' \
                        -f bam \
                        -l 1 \
                        -o {output.bam_cons} \
                        /dev/stdin \
                        2> {log.samb_log}
                """

rule sort_cons:
    input:
        bam_cons = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.bam"
    output:
        bam_cons_sort = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}.cons.sort.bam"
    log:
        sort_cons_log = config['resultsdir'] + "/logs/8_dd-cfDNA/{sample}.cons.sort.sambamba.log"
    params:
        temp_dir = config['resultsdir'] + "/results/tmp/"
    conda:
        "envs/twist_target.yaml"
    threads: sort_filt_threads
    shell:
        """
        sambamba sort \
            -t {threads} \
            --tmpdir {params.temp_dir} \
            -o {output.bam_cons_sort} \
            -l 1 \
            {input.bam_cons} \
        2> {log.sort_cons_log}
        """

