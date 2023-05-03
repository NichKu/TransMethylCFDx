rule demux_mv_UMI:
    input:
        samplesheet = config["runfolder"] + "/SampleSheet.csv"
    output:
        config["out_fastq"] + "/Stats/Stats.json"
    threads: demultiplex_threads
    params:
        runfolder = config["runfolder"],
        dir = config["out_fastq"]
    shell:
        """
        srslyumi {params.runfolder} {params.dir}
        """