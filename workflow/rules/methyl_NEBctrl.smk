"""
extract the controls lambda and pUC19
"""

rule extract_lambda_bam:
    input:
        bam_wCEREBIS = config['resultsdir'] + "/results/4_alignment/bam/{sample}_wCEREBIS_sorted.bam",
        bam_sorted_bai = config['resultsdir'] + "/results/4_alignment/bam/{sample}_wCEREBIS_sorted.bam.bai"
    output:
        bam_lambda = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda.bam"
    log:
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        samtools view -b -h -o {output.bam_lambda} {input.bam_wCEREBIS} "J02459.1"
        """

rule extract_pUC19_bam:
    input:
        bam_wCEREBIS = config['resultsdir'] + "/results/4_alignment/bam/{sample}_wCEREBIS_sorted.bam",
        bam_sorted_bai = config['resultsdir'] + "/results/4_alignment/bam/{sample}_wCEREBIS_sorted.bam.bai"
    output:
        bam_pUC19 = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19.bam"
    log:
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        samtools view -b -h -o {output.bam_pUC19} {input.bam_wCEREBIS} "M77789.2"
        """

"""
filter the bam files of the controls lambda and pUC19
"""

rule filter_bam_lambda:
    input:
        bam_lambda = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda.bam"
    output:
        bam_lambda_filt = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt.bam"
    log:
        filt_bam_log = config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}.filt_sambamba.log"
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        sambamba view \
            -h \
            --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam \
            -l 1 \
            -o {output.bam_lambda_filt} \
            {input.bam_lambda} \
        2> {log.filt_bam_log}
        """

rule filter_bam_pUC19:
    input:
        bam_pUC19 = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19.bam"
    output:
        bam_pUC19_filt = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt.bam"
    log:
        filt_bam_log = config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}.filt_sambamba.log"
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        sambamba view \
            -h \
            --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam \
            -l 1 \
            -o {output.bam_pUC19_filt} \
            {input.bam_pUC19} \
        2> {log.filt_bam_log}
        """

"""
move the UMI tag, index, sort and deduplicate the bam files of the controls lambda and pUC19
"""

rule move_umi_to_rx_tag_lambda:
        input:
                bam_lambda_filt = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt.bam"
        output:
                bam_sorted_UMIrx = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt.rx.bam"
        log:
                bamtag_log = config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}.rxbamtag.log"
        conda:
                "../envs/twist_target.yaml"
        shell:
                """
                srslyumi-bamtag \
                        --binary \
                        -o {output.bam_sorted_UMIrx} \
                        {input.bam_lambda_filt} \
                2> {log.bamtag_log}
                """

rule sort_bam_lambda:
        input:
                bam_lambda_filt = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt.rx.bam"
        output:
                bam_lambda_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_sort.rx.bam"
        log:
                config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}.sort.log"
        conda:
                "../envs/twist_target.yaml"
        threads: dedup_threads
        shell:
                """
                samtools sort \
                        -@ {threads} \
                        -o {output.bam_lambda_filt_sort} \
                        {input.bam_lambda_filt} \
                2> {log}
                """

rule index_umi_bam_lambda:
        input:
                bam_lambda_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_sort.rx.bam"
        output:
                config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_sort.rx.bam.bai"
        log:
                config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}.index_umibam.log"
        threads: index_threads
        shell:
                """
                samtools index -@ {threads} {input.bam_lambda_filt_sort}
                """

rule UMI_dedup_lambda:
        input:
                bam_lambda_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_sort.rx.bam",
                bam_lambda_filt_sort_bai = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_sort.rx.bam.bai"
        output:
                config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_edit_distance.tsv",
                config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_per_umi_per_position.tsv",
                config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_per_umi.tsv",
                bam_dedup = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_dedup_umitools.bam"
        log:
                umitoolsdedup_log = config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}.dedup_umitools.log"
        params:
                grouping_method = config['umi_tools_grouping_method'],
                output_prefix = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}"
        conda:
                "../envs/twist_target.yaml"
        shell:
                """
                umi_tools dedup \
                        --stdin={input.bam_lambda_filt_sort} \
                        --output-stats={params.output_prefix} \
                        --log={log.umitoolsdedup_log} \
                        --extract-umi-method=tag \
                        --umi-tag=RX \
                        --method {params.grouping_method} \
                        --paired \
                > {output.bam_dedup}
                """

rule move_umi_to_rx_tag_pUC19:
        input:
                bam_pUC19_filt = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt.bam"
        output:
                bam_sorted_UMIrx = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt.rx.bam"
        log:
                bamtag_log = config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}.rxbamtag.log"
        conda:
                "../envs/twist_target.yaml"
        shell:
                """
                srslyumi-bamtag \
                        --binary \
                        -o {output.bam_sorted_UMIrx} \
                        {input.bam_pUC19_filt} \
                2> {log.bamtag_log}
                """

rule sort_bam_pUC19:
        input:
                bam_pUC19_filt = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt.rx.bam"
        output:
                bam_pUC19_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_sort.rx.bam"
        log:
                config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}.sort.log"
        conda:
                "../envs/twist_target.yaml"
        threads: dedup_threads
        shell:
                """
                samtools sort \
                        -@ {threads} \
                        -o {output.bam_pUC19_filt_sort} \
                        {input.bam_pUC19_filt} \
                2> {log}
                """

rule index_umi_bam_pUC19:
        input:
                bam_pUC19_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_sort.rx.bam"
        output:
                config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_sort.rx.bam.bai"
        log:
                config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}.index_umibam.log"
        threads: index_threads
        shell:
                """
                samtools index -@ {threads} {input.bam_pUC19_filt_sort}
                """

rule UMI_dedup_pUC19:
        input:
                bam_pUC19_filt_sort = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_sort.rx.bam",
                bam_pUC19_filt_sort_bai = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_sort.rx.bam.bai"
        output:
                config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_edit_distance.tsv",
                config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_per_umi_per_position.tsv",
                config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_per_umi.tsv",
                bam_dedup = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_dedup_umitools.bam"
        log:
                umitoolsdedup_log = config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}.dedup_umitools.log"
        params:
                grouping_method = config['umi_tools_grouping_method'],
                output_prefix = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}"
        conda:
                "../envs/twist_target.yaml"
        shell:
                """
                umi_tools dedup \
                        --stdin={input.bam_pUC19_filt_sort} \
                        --output-stats={params.output_prefix} \
                        --log={log.umitoolsdedup_log} \
                        --extract-umi-method=tag \
                        --umi-tag=RX \
                        --method {params.grouping_method} \
                        --paired \
                > {output.bam_dedup}
                """

rule qc_methylation_lambda:
    input:
        bam_dedup = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}_lambda_filt_dedup_umitools.bam"
    output:
        config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}/{sample}_CpHRetentionByReadPos.txt",
        config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}/{sample}_totalReadConversionRate.txt",
        config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}/{sample}_CpGRetentionByReadPos.txt",
        bsstrand = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}/{sample}_strand_qc.txt"
    params:
        reference = config['reference_lambda'],
        output_dir = config['resultsdir'] + "/results/9_methyl_ctrl/lambda/{sample}",
        assets = config['biscuit_assets_lambda_dir']
    log:
        qc_sh = config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}_lambda.biscuit_qc.log",
        biscuit_log = config['resultsdir'] + "/logs/9_methyl_ctrl/lambda/{sample}_lambda.biscuit_bsstrand.log"
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        QC.sh \
            --outdir {params.output_dir} \
            {params.assets} \
            {params.reference} \
            {wildcards.sample} \
            {input.bam_dedup} \
        2> {log.qc_sh} || true
        biscuit bsstrand \
            {params.reference} \
            {input.bam_dedup} \
        1> {output.bsstrand}
        2> {log.biscuit_log}
        """

rule qc_methylation_puc19:
    input:
        bam_dedup = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}_pUC19_filt_dedup_umitools.bam"
    output:
        config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}/{sample}_CpHRetentionByReadPos.txt",
        config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}/{sample}_totalReadConversionRate.txt",
        config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}/{sample}_CpGRetentionByReadPos.txt",
        bsstrand = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}/{sample}_strand_qc.txt"
    params:
        reference = config['reference_pUC19'],
        output_dir = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19/{sample}",
        assets = config['biscuit_assets_pUC19_dir']
    log:
        qc_sh = config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}_pUC19.biscuit_qc.log",
        biscuit_log = config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19/{sample}_pUC19.biscuit_bsstrand.log"
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        QC.sh \
            --outdir {params.output_dir} \
            {params.assets} \
            {params.reference} \
            {wildcards.sample} \
            {input.bam_dedup} \
        2> {log.qc_sh} || true
        biscuit bsstrand \
            {params.reference} \
            {input.bam_dedup} \
        1> {output.bsstrand}
        2> {log.biscuit_log}
        """