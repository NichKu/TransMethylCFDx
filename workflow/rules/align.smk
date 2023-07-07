rule align:
    input:
        reads_1_trimmed = config['resultsdir'] + "/results/2_trimmed/{sample}_R1_val_1.fq.gz",
        reads_2_trimmed = config['resultsdir'] + "/results/2_trimmed/{sample}_R2_val_2.fq.gz"
    output:
        bam_w_ctrls = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls.bam"
    log:
        align_log = config['resultsdir'] + "/logs/4_alignment/{sample}.bwameth.log",
        samb_log = config['resultsdir'] + "/logs/4_alignment/{sample}.toBam_sambamba.log"
    params:
        reference = config['reference_w_ctrl']
    conda:
        "../envs/twist_target.yaml"
    threads: align_threads
    shell:
        """
        bwameth.py \
            --reference {params.reference} \
            -t {threads} \
            {input.reads_1_trimmed} {input.reads_2_trimmed} \
        2> {log.align_log} | \
        sambamba view \
            -h \
            -t {threads} \
            --sam-input \
            -f bam \
            -l 1 \
            -o {output.bam_w_ctrls} \
            /dev/stdin \
        2> {log.samb_log}
        """

rule sort_bam_w_ctrls:
    input:
        bam_w_ctrls = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls.bam"
    output:
        bam_w_ctrls_sort = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls_sorted.bam"
    log:
        sort_unfilt_log = config['resultsdir'] + "/logs/4_alignment/{sample}.sort_w_ctrls_sambamba.log"
    params:
        temp_dir = config['resultsdir'] + "/results/tmp/"
    conda:
        "../envs/twist_target.yaml"
    threads: sort_filt_threads
    shell:
        """
        sambamba sort \
            -t {threads} \
            --tmpdir {params.temp_dir} \
            -o {output.bam_w_ctrls_sort} \
            -l 1 \
            {input.bam_w_ctrls} \
        2> {log.sort_unfilt_log}
        """

rule index_bam_w_ctrls:
    input:
        bam_w_ctrls_sort = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls_sorted.bam"
    output:
        bam_w_ctrls_sort_bai = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls_sorted.bam.bai"
    log:
        log_index = config['resultsdir'] + "/logs/4_alignment/{sample}.index_bam_w_ctrls.log"
    conda:
        "../envs/twist_target.yaml"
    threads: index_threads
    shell:
        """
        samtools index \
            -@ {threads} \
            -o {output.bam_w_ctrls_sort_bai} \
            {input.bam_w_ctrls_sort} \
        2> {log.log_index}
        """

rule remove_CEREBIS:
    input:
        bam_w_ctrls_sort = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls_sorted.bam",
        bam_w_ctrls_sort_bai = config['resultsdir'] + "/results/4_alignment/bam/{sample}_w_ctrls_sorted.bam.bai"
    output:
        bam = config['resultsdir'] + "/results/4_alignment/bam/{sample}_woCEREBIS.bam"
    log:
        log_remove = config['resultsdir'] + "/logs/4_alignment/{sample}.removeCEREBIS.log"
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        samtools view \
            -b \
            -h \
            -o {output.bam} \
            {input.bam_w_ctrls_sort} \
            chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
        2> {log.log_remove}
        """

rule sort_unfilt_woCEREBIS:
    input:
        bam_unfilt = config['resultsdir'] + "/results/4_alignment/bam/{sample}_woCEREBIS.bam"
    output:
        bam_unfilt_sort = config['resultsdir'] + "/results/4_alignment/bam/{sample}_unfilt_sorted.bam"
    log:
        sort_unfilt_log = config['resultsdir'] + "/logs/4_alignment/{sample}.sort_unfilt_sambamba.log"
    params:
        temp_dir = config['resultsdir'] + "/results/tmp/"
    conda:
        "../envs/twist_target.yaml"
    threads: sort_filt_threads
    shell:
        """
        sambamba sort \
            -t {threads} \
            --tmpdir {params.temp_dir} \
            -o {output.bam_unfilt_sort} \
            -l 1 \
            {input.bam_unfilt} \
        2> {log.sort_unfilt_log}
        """


rule filter_bam_woCEREBIS:
    input:
        bam_unfilt = config['resultsdir'] + "/results/4_alignment/bam/{sample}_woCEREBIS.bam"
    output:
        bam_filt =  config['resultsdir'] + "/results/4_alignment/bam_filt/{sample}_filt.bam"
    log:
        filt_bam_log = config['resultsdir'] + "/logs/4_alignment/{sample}.filt_sambamba.log",
        sort_filt_log = config['resultsdir'] + "/logs/4_alignment/{sample}.sort_filt_sambamba.log"
    params:
        temp_dir = config['resultsdir'] + "/results/tmp/"
    conda:
        "../envs/twist_target.yaml"
    threads: sort_filt_threads
    shell:
        """
        sambamba view \
            -h \
            -t {threads} \
            --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam \
            -l 1 \
            -o /dev/stdout \
            {input.bam_unfilt} \
        2> {log.filt_bam_log} | \
        sambamba sort \
            -t {threads} \
            -o {output.bam_filt} \
            --tmpdir {params.temp_dir} \
            -l 1 \
            /dev/stdin \
        2> {log.sort_filt_log}
        """

rule index_filt_bam:
    input:
        bam_sorted = config['resultsdir'] + "/results/4_alignment/bam_filt/{sample}_filt.bam"
    output:
        bam_sorted_bai = config['resultsdir'] + "/results/4_alignment/bam_filt/{sample}_filt.bam.bai"
    log:
        config['resultsdir'] + "/logs/4_alignment/{sample}.index_filt_bam.log"
    conda:
        "../envs/twist_target.yaml"
    threads: index_threads
    shell:
        """
        samtools index -@ {threads} -o {output.bam_sorted_bai} {input.bam_sorted}
        """