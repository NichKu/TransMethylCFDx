rule fastqc_trimmed:
    input:
        reads_1_trimmed = config['resultsdir'] + "/results/2_trimmed/{sample}_R1_val_1.fq.gz",
        reads_2_trimmed = config['resultsdir'] + "/results/2_trimmed/{sample}_R2_val_2.fq.gz"
    output:
        html_1 = config['resultsdir'] + "/results/3_fastqc_trimmed/{sample}_R1_val_1_fastqc.html",
        zip_1 = config['resultsdir'] + "/results/3_fastqc_trimmed/{sample}_R1_val_1_fastqc.zip",
        html_2 = config['resultsdir'] + "/results/3_fastqc_trimmed/{sample}_R2_val_2_fastqc.html",
        zip_2 = config['resultsdir'] + "/results/3_fastqc_trimmed/{sample}_R2_val_2_fastqc.zip"
    log:
        fastqc_trimmed_log = config['resultsdir'] + "/logs/3_fastqc_trimmed/{sample}.fastqc_trimmed.log"
    conda:
        "../envs/twist_target.yaml"
    params:
        output_dir = config['resultsdir'] + "/results/3_fastqc_trimmed/"
    threads: fastqc_threads
    shell:
        """
        fastqc --threads {threads} -o {params.output_dir} {input.reads_1_trimmed} {input.reads_2_trimmed} 2> {log.fastqc_trimmed_log}
        """
