rule multiqc:
    input:
        expand(config['resultsdir'] + "/results/1_fastqc_init/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/3_fastqc_trimmed/{sample}_R1_val_1_fastqc.zip", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/5_dedup/{sample}_filt.umitools.tsv", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/5_dedup/{sample}{ending}", 
            ending=["_edit_distance.tsv", "_per_umi_per_position.tsv", "_per_umi.tsv"], sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/meth_metrics/{sample}.meth_hs_metrics.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.meth_per_target_coverage.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/snp_metrics/{sample}.snp_hs_metrics.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.snp_per_target_coverage.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.bam.flagstat", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.lc_extrap.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/{sample}/{sample}.idxstats.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/{sample}/qualimapReport.html", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_CpGRetentionByReadPos.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_strand_qc.txt", sample=SAMPLES)
    output:
        multiqc_report = config['resultsdir'] + "/results/10_multiqc/multiqc_report.html",
        multiqc_report_cerebis = config['resultsdir'] + "/results/10_multiqc/multiqc_report_CEREBIS.html"
    log:
        multiqc_log_res = config['resultsdir'] + "/logs/10_multiqc/multiqc_report.results.log",
        multiqc_log_cerebis = config['resultsdir'] + "/logs/10_multiqc/multiqc_report.cerebis.log"
    params:
        scan_directory = config['resultsdir'],
        scan_directory_exlusion_1 = config['resultsdir'] + "/results/9_methyl_ctrl",
        scan_directory_exlusion_2 = config['resultsdir'] + "/logs/9_methyl_ctrl",
        scan_directory_cerebis_res = config['resultsdir'] + "/results/9_methyl_ctrl",
        scan_directory_cerebis_log = config['resultsdir'] + "/logs/9_methyl_ctrl",
        out_directory = config['resultsdir'] + "/results/10_multiqc/"
    conda:
        "envs/twist_target.yaml"
    threads: finalQC_threads
    shell:
        """
        multiqc {params.scan_directory} -f -d -dd 1 --ignore {params.scan_directory_exlusion_1} --ignore {params.scan_directory_exlusion_2} --outdir {params.out_directory} 2> {log.multiqc_log_res}
        #multiqc {params.scan_directory_cerebis_res} {params.scan_directory_cerebis_log} -f --outdir {params.out_directory} --filename multiqc_report_CEREBIS.html 2> {log.multiqc_log_cerebis}
        """