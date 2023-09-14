rule trim_reads:
#   Description: trims the provided raw reads
#   todo: switch trimmomatic for fastp?
#   todo: implement autodetection of compositional bias trimming?
#   todo: do we provide the adapter list? switching to fastp would provide automatic adapter identification
    input:
        r1 = config['R1']
        r2 = config['R2']
        adapters = config['adapters']
    output:
        r1 = "trimmed_reads/" + config['R1'].split("/")[-1]
        r1_unpaired = "trimmed_reads/unpaired." + config['R1'].split("/")[-1]
        r2 = "trimmed_reads/" + config['R2'].split("/")[-1]
        r2_unpaired = "trimmed_reads/unpaired." + config['R2'].split("/")[-1]
    threads: config['threads']
    shell:
        """
        mkdir -p trimmed_reads
        trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} ILLUMINACLIP:{input.adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15"
        """

rule assemble_transcriptome:
#   Description: Assembles a transcriptome if it is not provided. 
#   In this case sequencing reads MUST be provided in the config.
    input:
        r1 = rules.trim_reads.output.r1
        r2 = rules.trim_reads.output.r2
    

rule cdhit_clustering:
#   Description: Clusters transcriptome sequences using cd-hit-est.
    input: 
        transcriptome = config['transcriptome']
    output:
        clustered_transcriptome = "clustered.{input.transcriptome}"
    params:
        threshold = config['clustering_threshold']
        memory = str(int(config['memory'])*1000)
    threads: config['threads']
    shell:
        """
        cd-hit-est -i {input.transcriptome} -o {output.clustered_transcriptome} -c {params.threshold} -M {params.memory} -T {threads} 
        """

