rule fastqc_init:
    input:
        reads_1=config["fastqdir"] + "/{sample}_R1.fastq.gz",
        reads_2=config["fastqdir"] + "/{sample}_R2.fastq.gz",
    output:
        html_1=config["resultsdir"] + "/results/1_fastqc_init/{sample}_R1_fastqc.html",
        zip_1=config["resultsdir"] + "/results/1_fastqc_init/{sample}_R1_fastqc.zip",
        html_2=config["resultsdir"] + "/results/1_fastqc_init/{sample}_R2_fastqc.html",
        zip_2=config["resultsdir"] + "/results/1_fastqc_init/{sample}_R2_fastqc.zip",
    log:
        fastqc_init_log=config["resultsdir"]
        + "/logs/1_fastqc_init/{sample}.fastqc_init.log",
    conda:
        "../envs/twist_target.yaml"
    params:
        output_dir=config["resultsdir"] + "/results/1_fastqc_init/",
    threads: fastqc_threads
    shell:
        """
        fastqc --threads {threads} -o {params.output_dir} {input.reads_1} {input.reads_2} 2> {log.fastqc_init_log}
        """
