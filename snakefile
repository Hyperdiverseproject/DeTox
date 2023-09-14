
rule cdhit_clustering:
#   Description: Clusters transcriptome sequences using cd-hit-est.
    input: 
        transcriptome = config['transcriptome']
    output:
        clustered_transcriptome = "{basename}.clustered.fasta"
    params:
        threshold = config['clustering_threshold']
        memory = str(int(config['memory'])*1000)
        threads = config['threads']
    threads: config['threads']
    shell:
        """
        cd-hit-est -i {input.transcriptome} -o {output.clustered_transcriptome} -c {params.threshold} -M {params.memory} -T {params.threads} 
        """

