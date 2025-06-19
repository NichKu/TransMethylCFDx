rule trim_adapter:
    input:
        reads_1=config["fastqdir"] + "/{sample}_R1.fastq.gz",
        reads_2=config["fastqdir"] + "/{sample}_R2.fastq.gz",
    output:
        reads_1_trimmed=config["resultsdir"]
        + "/results/2_trimmed/{sample}_R1_val_1.fq.gz",
        reads_2_trimmed=config["resultsdir"]
        + "/results/2_trimmed/{sample}_R2_val_2.fq.gz",
    log:
        trim_log=config["resultsdir"] + "/logs/2_trimmed/{sample}.trim.log",
    params:
        output_dir=config["resultsdir"] + "/results/2_trimmed/",
        #trim_galore_clip_r1 = config['trim_galore_clip_r1'],
        #trim_galore_clip_r2 = config['trim_galore_clip_r2']
    conda:
        "../envs/twist_target.yaml"
    threads: trim_threads
    shell:
        """
        trim_galore \
            --gzip \
            --cores {threads} \
            --output_dir {params.output_dir} \
            --2colour 20 \
            --paired {input.reads_1} {input.reads_2} \
        2> {log.trim_log}
        """


# --clip_R1 {params.trim_galore_clip_r1} \
# --clip_R2 {params.trim_galore_clip_r2} \
