rule qc_alignment_dup:
    input:
        bam_unfilt_sort = config['resultsdir'] + "/results/4_alignment/bam/{sample}_unfilt_sorted.bam",
        bam_filt =  config['resultsdir'] + "/results/4_alignment/bam_filt/{sample}_filt.bam"
    output:
        samtools_flagstat = config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.bam.flagstat",
        preseq_estimate = config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.lc_extrap.txt"
    log:
        samtools_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.samtools.log",
        preseq_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.preseq.log"
    conda:
        "envs/twist_target.yaml"
    threads: qc_align_threads
    shell:
        """
        samtools stats -@ {threads} {input.bam_unfilt_sort} 1> {output.samtools_flagstat} 2> {log.samtools_log}
        preseq lc_extrap -B -P -o {output.preseq_estimate} {input.bam_filt} 2> {log.preseq_log}
        """

rule qc_aligment_dedup:
    input:
        bam_dedup = config['resultsdir'] + "/results/5_dedup/{sample}_filt_dedup_umitools.bam",
        all_targets_intervallist = output_file_all_target_interval_list,
        meth_targets_intervallist = output_file_meth_target_interval_list,
        snp_targets_intervallist = output_file_snp_target_interval_list,
        all_bait_intervallist = output_file_all_bait_interval_list,
        meth_bait_intervallist = output_file_meth_bait_interval_list,
        snp_bait_intervallist = output_file_snp_bait_interval_list
    output:
        config['resultsdir'] + "/results/6_aligment_QC/{sample}/qualimapReport.html",
        all_hs_metrics = config['resultsdir'] + "/results/6_aligment_QC/all_metrics/{sample}.all_hs_metrics.txt",
        meth_hs_metrics = config['resultsdir'] + "/results/6_aligment_QC/meth_metrics/{sample}.meth_hs_metrics.txt",
        meth_per_targ_cov = config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.meth_per_target_coverage.txt",
        snp_hs_metrics = config['resultsdir'] + "/results/6_aligment_QC/snp_metrics/{sample}.snp_hs_metrics.txt",
        snp_per_targ_cov = config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.snp_per_target_coverage.txt",

        idxstats = config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.idxstats.txt"
    log:
        all_hsmetrics_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.all_hsmetrics.log",
        meth_hsmetrics_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.meth_hsmetrics.log",
        snp_hsmetrics_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.snp_hsmetrics.log",
        idxstats_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.idxstats.log",
        qualimap_log = config['resultsdir'] + "/logs/6_aligment_QC/{sample}.bamqc.log"
    params:
        qualimap_outdir = config['resultsdir'] + "/results/6_aligment_QC/{sample}",
        panel_region = config['targetregions_all'],
        ref = config["reference_w_ctrl"]
    conda:
        "envs/twist_target.yaml"
    resources:
        mem_mb=6000
    threads: qc_align_threads
    shell:
        """
        picard CollectHsMetrics \
            I={input.bam_dedup} \
            O={output.all_hs_metrics} \
            R={params.ref} \
            TARGET_INTERVALS={input.all_targets_intervallist} \
            BAIT_INTERVALS={input.all_bait_intervallist} \
            COVERAGE_CAP=10000 \
        2> {log.all_hsmetrics_log}

        picard CollectHsMetrics \
            I={input.bam_dedup} \
            O={output.meth_hs_metrics} \
            R={params.ref} \
            TARGET_INTERVALS={input.meth_targets_intervallist} \
            BAIT_INTERVALS={input.meth_bait_intervallist} \
            PER_TARGET_COVERAGE={output.meth_per_targ_cov} \
            COVERAGE_CAP=10000 \
        2> {log.meth_hsmetrics_log}

        picard CollectHsMetrics \
            I={input.bam_dedup} \
            O={output.snp_hs_metrics} \
            R={params.ref} \
            TARGET_INTERVALS={input.snp_targets_intervallist} \
            BAIT_INTERVALS={input.snp_bait_intervallist} \
            PER_TARGET_COVERAGE={output.snp_per_targ_cov} \
            COVERAGE_CAP=10000 \
        2> {log.snp_hsmetrics_log}
        
        samtools idxstats \
            -@ {threads} \
            {input.bam_dedup} \
        1> {output.idxstats} \
        2> {log.idxstats_log}

        qualimap bamqc --java-mem-size={resources.mem_mb}M \
            -gff {params.panel_region} \
            -ip \
            -nt {threads} \
            -os \
            -sdmode 1 \
            -bam {input.bam_dedup} \
            -outdir {params.qualimap_outdir} \
            -outformat HTML \
        2> {log.qualimap_log}
        """

rule qc_methylation:
    input:
        bam_dedup = config['resultsdir'] + "/results/5_dedup/{sample}_filt_dedup_umitools.bam"
    output:
        config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_totalReadConversionRate.txt",
        config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_CpHRetentionByReadPos.txt",
        config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_CpGRetentionByReadPos.txt",
        bsstrand = config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_strand_qc.txt"
    params:
        reference = config['reference_wo_ctrl'],
        output_dir = config['resultsdir'] + "/results/7_methylQC/{sample}",
        assets = config['biscuit_assets_hg19_dir']
    log:
        qc_sh = config['resultsdir'] + "/logs/7_methylQC/{sample}.biscuit_qc.log",
        biscuit_log = config['resultsdir'] + "/logs/7_methylQC/{sample}.biscuit_bsstrand.log"
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