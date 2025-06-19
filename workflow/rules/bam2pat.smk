rule bam2pat:
    input:
        bam_cons=config["resultsdir"]
        + "/results/5_dedup_consensus/{sample}.cons.sort.bam",
    output:
        pat_cons=config["resultsdir"] + "/results/9_pat_files/{sample}.cons.sort.pat.gz",
    params:
        temp_dir=config["resultsdir"] + "/results/tmp/",
        output_dir=config["resultsdir"] + "/results/",
    log:
        bam2pat_log=config["resultsdir"] + "/logs/9_pat_files/",
    conda:
        "../envs/twist_target.yaml"
    threads: bam2pat_threads
    shell:
        """
        wgbstools bam2pat \
            -@ {threads.bam2pat_threads} \
            --lbeta \
            -o {params.output_dir} \
            -T {params.temp_dir} \
            {input.bam_cons}
        2> {log.bam2pat_log}
        """
