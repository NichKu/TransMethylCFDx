rule extract_cerebis_bam:
    input:
        bam_wCEREBIS = config['resultsdir'] + "/results/4_alignment/bam/{sample}_wCEREBIS_sorted.bam",
        bam_sorted_bai = config['resultsdir'] + "/results/4_alignment/bam/{sample}_wCEREBIS_sorted.bam.bai"
    output:
        bam_CEREBIS = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS.bam"
    log:
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        samtools view -b -h -o {output.bam_CEREBIS} {input.bam_wCEREBIS} CEREBIS
        """

rule filter_bam_cerebis:
    input:
        bam_CEREBIS = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS.bam"
    output:
        bam_CEREBIS_filt = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt.bam"
    log:
        filt_bam_log = config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}_CEREBIS.filt_sambamba.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        sambamba view \
            -h \
            --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam \
            -l 1 \
            -o {output.bam_CEREBIS_filt} \
            {input.bam_CEREBIS} \
        2> {log.filt_bam_log}
        """

rule move_umi_to_rx_tag_cerebis:
        input:
                bam_CEREBIS_filt = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt.bam"
        output:
                bam_sorted_UMIrx = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt.rx.bam"
        log:
                bamtag_log = config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}.rxbamtag.log"
        conda:
                "envs/twist_target.yaml"
        shell:
                """
                srslyumi-bamtag \
                        --binary \
                        -o {output.bam_sorted_UMIrx} \
                        {input.bam_CEREBIS_filt} \
                2> {log.bamtag_log}
                """

rule sort_bam_cerebis:
        input:
                bam_CEREBIS_filt = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt.rx.bam"
        output:
                bam_CEREBIS_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_sort.rx.bam"
        log:
                config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}.sort.log"
        conda:
                "envs/twist_target.yaml"
        threads: dedup_threads
        shell:
                """
                samtools sort \
                        -@ {threads} \
                        -o {output.bam_CEREBIS_filt_sort} \
                        {input.bam_CEREBIS_filt} \
                2> {log}
                """

rule index_umi_bam_cerebis:
        input:
                bam_CEREBIS_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_sort.rx.bam"
        output:
                config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_sort.rx.bam.bai"
        log:
                config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}.index_umibam.log"
        threads: index_threads
        shell:
                """
                samtools index -@ {threads} {input.bam_CEREBIS_filt_sort}
                """

rule UMI_dedup_cerebis:
        input:
                bam_CEREBIS_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_sort.rx.bam",
                bam_CEREBIS_filt_sort_bai = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_sort.rx.bam.bai"
        output:
                config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_edit_distance.tsv",
                config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_per_umi_per_position.tsv",
                config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_per_umi.tsv",
                bam_dedup = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_dedup_umitools.bam"
        log:
                umitoolsdedup_log = config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}.dedup_umitools.log"
        params:
                grouping_method = config['umi_tools_grouping_method'],
                output_prefix = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}"
        conda:
                "envs/twist_target.yaml"
        shell:
                """
                umi_tools dedup \
                        --stdin={input.bam_CEREBIS_filt_sort} \
                        --output-stats={params.output_prefix} \
                        --log={log.umitoolsdedup_log} \
                        --extract-umi-method=tag \
                        --umi-tag=RX \
                        --method {params.grouping_method} \
                        --paired \
                > {output.bam_dedup}
                """

rule qc_methylation_cerebis:
    input:
        bam_dedup = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}_CEREBIS_filt_dedup_umitools.bam"
    output:
        config['resultsdir'] + "/results/9_methyl_ctrl/{sample}/{sample}_CpHRetentionByReadPos.txt",
        config['resultsdir'] + "/results/9_methyl_ctrl/{sample}/{sample}_totalReadConversionRate.txt",
        config['resultsdir'] + "/results/9_methyl_ctrl/{sample}/{sample}_CpGRetentionByReadPos.txt",
        bsstrand = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}/{sample}_strand_qc.txt"
    params:
        reference = config['reference_ctrl'],
        output_dir = config['resultsdir'] + "/results/9_methyl_ctrl/{sample}",
        assets = config['biscuit_assets_cerebis_dir']
    log:
        qc_sh = config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}_CEREBIS.biscuit_qc.log",
        biscuit_log = config['resultsdir'] + "/logs/9_methyl_ctrl/{sample}_CEREBIS.biscuit_bsstrand.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        QC.sh \
            --outdir {params.output_dir} \
            {params.assets} \
            {params.reference} \
            {wildcards.sample} \
            {input.bam_dedup} \        
        2> {log.qc_sh}
        biscuit bsstrand \
            {params.reference} \
            {input.bam_dedup} \
        1> {output.bsstrand}
        2> {log.biscuit_log}
        """