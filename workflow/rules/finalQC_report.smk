rule multiqc:
    input:
        expand(config['resultsdir'] + "/results/1_fastqc_init/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/3_fastqc_trimmed/{sample}_R1_val_1_fastqc.zip", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/5_dedup_consensus/{sample}_filt.umitools.tsv", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/5_dedup_consensus/{sample}{ending}", 
            ending=["_edit_distance.tsv", "_per_umi_per_position.tsv", "_per_umi.tsv"], sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/flagstat/{sample}.bam.flagstat", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/preseq/{sample}.lc_extrap.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/all_metrics/{sample}.all_hs_metrics.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/meth_metrics/{sample}.meth_hs_metrics.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/meth_per_target_metrics/{sample}.meth_per_target_coverage.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/snp_metrics/{sample}.snp_hs_metrics.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/snp_per_target_metrics/{sample}.snp_per_target_coverage.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/gcbias_metrics/{sample}.{ext}", sample=SAMPLES, 
            ext=["bam_gc_metrics.txt", "bam_gcbias_summary.txt", "bam_gcbias.pdf"]),
        expand(config['resultsdir'] + "/results/6_aligment_QC/idxstats/{sample}.idxstats.txt", sample=SAMPLES),
        expand(config['resultsdir'] + "/results/6_aligment_QC/qualimap/{sample}/{sample}.qualimapReport.pdf", sample=SAMPLES)
        #expand(config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_CpGRetentionByReadPos.txt", sample=SAMPLES),
        #expand(config['resultsdir'] + "/results/7_methylQC/{sample}/{sample}_strand_qc.txt", sample=SAMPLES)
    output:
        multiqc_report = config['resultsdir'] + "/results/10_multiqc/multiqc_report_" + config["project_run_name"] + ".html",
        multiqc_report_cerebis = config['resultsdir'] + "/results/10_multiqc/multiqc_report_" + config["project_run_name"] + "_CEREBIS.html",
        #multiqc_report_lambda = config['resultsdir'] + "/results/10_multiqc/multiqc_report_" + config["project_run_name"] + "_lambda.html",
        #multiqc_report_pUC19 = config['resultsdir'] + "/results/10_multiqc/multiqc_report_" + config["project_run_name"] + "_pUC19.html"
    log:
        multiqc_log_res = config['resultsdir'] + "/logs/10_multiqc/multiqc_report.results.log",
        multiqc_log_cerebis = config['resultsdir'] + "/logs/10_multiqc/multiqc_report.cerebis.log"
        #multiqc_log_lambda = config['resultsdir'] + "/logs/10_multiqc/multiqc_report.lambda.log",
        #multiqc_log_pUC19 = config['resultsdir'] + "/logs/10_multiqc/multiqc_report.pUC19.log"
    params:
        scan_directory = config['resultsdir'],
        scan_directory_exlusion_1 = config['resultsdir'] + "/results/9_methyl_ctrl",
        scan_directory_exlusion_2 = config['resultsdir'] + "/logs/9_methyl_ctrl",
        scan_directory_cerebis_res = config['resultsdir'] + "/results/9_methyl_ctrl/CEREBIS",
        scan_directory_cerebis_log = config['resultsdir'] + "/logs/9_methyl_ctrl/CEREBIS",
        scan_directory_lambda_res = config['resultsdir'] + "/results/9_methyl_ctrl/lambda",
        scan_directory_lambda_log = config['resultsdir'] + "/logs/9_methyl_ctrl/lambda",
        scan_directory_pUC19_res = config['resultsdir'] + "/results/9_methyl_ctrl/pUC19",
        scan_directory_pUC19_log = config['resultsdir'] + "/logs/9_methyl_ctrl/pUC19",
        out_directory = config['resultsdir'] + "/results/10_multiqc/",
        multiqc_config = config['multiqc_config'],
        project_name = config["project_run_name"]
    conda:
        "../envs/twist_target.yaml"
    threads: finalQC_threads
    shell:
        """
        pwd
        multiqc {params.scan_directory} -f -d -dd 1 --ignore {params.scan_directory_exlusion_1} --ignore {params.scan_directory_exlusion_2} --outdir {params.out_directory} --filename multiqc_report_{params.project_name} 2> {log.multiqc_log_res}
        multiqc {params.scan_directory_cerebis_res} {params.scan_directory_cerebis_log} -f --outdir {params.out_directory} --filename multiqc_report_{params.project_name}_CEREBIS.html 2> {log.multiqc_log_cerebis}
        """
#-c {params.multiqc_config} # creates "The 'picard' MultiQC module broke..." error
#multiqc {params.scan_directory_lambda_res} {params.scan_directory_lambda_log} -f --outdir {params.out_directory} --filename multiqc_report_{params.project_name}_lambda.html 2> {log.multiqc_log_lambda}
#multiqc {params.scan_directory_pUC19_res} {params.scan_directory_pUC19_log} -f --outdir {params.out_directory} --filename multiqc_report_{params.project_name}_pUC19.html 2> {log.multiqc_log_pUC19}