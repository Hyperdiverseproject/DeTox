rule trim_reads:
#   Description: trims the provided raw reads
#   todo: switch trimmomatic for fastp?
#   todo: implement autodetection of compositional bias trimming?
#   todo: do we provide the adapter list? switching to fastp would provide automatic adapter identification
    input:
        r1 = config['R1'],
        r2 = config['R2'],
        adapters = config['adapters'],
    output:
        r1 = "trimmed_reads/" + config['R1'].split("/")[-1],
        r1_unpaired = "trimmed_reads/unpaired." + config['R1'].split("/")[-1],
        r2 = "trimmed_reads/" + config['R2'].split("/")[-1],
        r2_unpaired = "trimmed_reads/unpaired." + config['R2'].split("/")[-1],
    threads: config['threads']
    shell:
        """
        mkdir -p trimmed_reads
        trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} ILLUMINACLIP:{input.adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15"
        """

rule assemble_transcriptome:
#   Description: Assembles a transcriptome if it is not provided. Uses Trinity 
#   In this case sequencing reads MUST be provided in the config.
    input:
        r1 = rules.trim_reads.output.r1,
        r2 = rules.trim_reads.output.r2,
    output:
        assembly = "trinity_out_dir/Trinity.fasta"
    params:
        memory = config['memory'] + "G"
    threads: config['threads']
    shell:
        """
        Trinity --seqType fq --left {rules.trim_reads.output.r1} --right {rules.trim_reads.output.r2} --CPU {threads} --max_memory {params.memory} --KMER_SIZE 31
        """
    

rule cdhit_clustering:
#   Description: Clusters transcriptome sequences using cd-hit-est.
    input: 
        transcriptome = config['transcriptome'] if transcriptome in config else rules.assemble_transcriptome.output.assembly
    output:
        clustered_transcriptome = config['basename'] + ".clustered.fasta",
        basename = config['basename'] + ".clustered",
    params:
        threshold = config['clustering_threshold'],
        memory = str(int(config['memory'])*1000),
    threads: config['threads']
    shell:
        """
        cd-hit-est -i {input.transcriptome} -o {output.basename} -c {params.threshold} -M {params.memory} -T {threads} 
        """

rule build_contaminants_database:
#   Description: builds blast database for the removal of contaminants   
#   todo: make this optional like in the original code
    input:
        fasta_db = config['contaminants']
    output:
        blast_db = config['contaminants'] + ".nin"
    shell:
        """
        makeblastdb -dbtype nucl -in {input.fasta_db}
        """

rule blast_on_contaminants:
#   Description: performs the actual blast of the contigs against the contaminants database
    input:
        blast_db = rules.build_contaminants_database.output.blast_db,
        contigs = rules.cdhit_clustering.output.clustered_transcriptome,
    output:
        blast_result = config['basename'] + ".blastsnuc.out"
    params:
        evalue = config['contamination_evalue']
    threads: config['threads']
    shell:
        """
        blastn -db {input.blast_db} -query {input.contigs} -out {output.blast_result} -outfmt 6 -evalue {params.evalue} -max_target_seqs 1 -numthreads {threads}
        """

rule filter_contaminants:
#   Description: performs the actual filtering
    input: 
        blast_result = rules.blast_on_contaminants.output.blast_result,
        contigs = rules.cdhit_clustering.output.clustered_transcriptome,
    output:
        filtered_contigs = config['basename'] + ".filtered.fasta"
    run:
        from Bio import SeqIO
        records = []
        infile = open(input.blast_result, 'r')
        for line in infile:
            line = line.rstrip()
            if line[0] != '#':
                blast = line.split()                                        
                records.append(blast[0]) # we recover the ID of the significan hits
        infile.close()
        recordIter = SeqIO.parse(open(input.contigs), "fasta")
        with open(output.filtered_contigs, "w") as handle:
            for rec in recordIter:
                if rec.id not in listIDConta:
                    SeqIO.write(rec, handle, "fasta")
