def get_signalp_splits(wildcards):
    import os
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    return expand('./split_files/{i}.fasta',
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.fasta')).i)



#configfile: "config.yaml"



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
        memory = str(config['memory']) + "G"
    threads: config['threads']
    shell:
        """
        Trinity --seqType fq --left {rules.trim_reads.output.r1} --right {rules.trim_reads.output.r2} --CPU {threads} --max_memory {params.memory} --KMER_SIZE 31
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
        blast_db_alias =  config['contaminants'],
        contigs = config['transcriptome'] if 'transcriptome' in config else rules.assemble_transcriptome.output.assembly,
    output:
        blast_result = config['basename'] + ".blastsnuc.out"
    params:
        evalue = config['contamination_evalue']
    threads: config['threads']
    shell:
        """
        blastn -db {input.blast_db_alias} -query {input.contigs} -out {output.blast_result} -outfmt 6 -evalue {params.evalue} -max_target_seqs 1 -num_threads {threads}
        """

rule filter_contaminants:
#   Description: performs the actual filtering
    input: 
        blast_result = rules.blast_on_contaminants.output.blast_result,
        contigs = config['transcriptome'] if 'transcriptome' in config else rules.assemble_transcriptome.output.assembly,
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
                if rec.id not in records:
                    SeqIO.write(rec, handle, "fasta")

rule detect_orfs:
#   Description: finds complete orfs within the input nucleotide sequences. 
#   i'm testing this with orfipy instead of orffinder to leverage multithreading
    input:
        nucleotide_sequences = rules.filter_contaminants.output.filtered_contigs
    output:
        aa_sequences = config['basename'] + ".faa"
    params:
        minlen = "66" if "minlen" not in config else config['minlen'],
        maxlen = "30000000" if "maxlen" not in config else config['maxlen']
    threads: config['threads']
    shell:
        """
        orfipy --procs {threads} --start ATG --pep {output.aa_sequences} --min {params.minlen} --max {params.maxlen} {input.nucleotide_sequences} --outdir .
        """

rule cluster_peptides:
#   Description: runs cd-hit on predicted peptide to remove excess redundancy
    input:
        aa_sequences = rules.detect_orfs.output.aa_sequences
    output:
        filtered_aa_sequences = config['basename'] + ".filtered.faa" #todo might want to change the basename to a param also in the other cd-hit rule if we decide on keeping it
    params:
        threshold = config['clustering_threshold'], #todo : we might want to use separate thresholds if we're going to run cd-hit on transcripts and peptides
        memory = str(int(config['memory'])*1000),
        basename = config['basename'] + ".filtered"
    threads: config['threads']
    shell:
        """
        cd-hit -i {input.aa_sequences} -o {output.filtered_aa_sequences} -c {params.threshold} -M {params.memory} -T {threads} 
        """

rule trim_peptides:
#   Description: this rule trims all peptides to only the first 50 aminoacids, as they are the only useful part for signalp. This step improves load time.
    input:
        aa_sequences = rules.cluster_peptides.output.filtered_aa_sequences
    output:
        trimmed_sequences = config['basename'] + ".trimmed.faa",
    threads: 
        config['threads']
    run:
        from Bio import SeqIO
        import subprocess
        with open(output.trimmed_sequences, "w") as outfile:
            for seq in SeqIO.parse(input.aa_sequences, "fasta"):
                outfile.write(f">{seq.id}\n{seq.seq[:31]}\n")


checkpoint split_fasta:
    input:
        fasta_file = rules.trim_peptides.output.trimmed_sequences
    output:
        split_dir = directory("split_files"),
        marker = "split_fasta.done"
    params:
        chunk_size = 9000 # using 9000 instead of 50000 for usability in normal desktop/laptop pcs. May be user defined.
    run:
        from Bio import SeqIO
        import os
        def batch_iterator(iterator, batch_size):
            """Returns lists of length batch_size.

            This is a generator function, and it returns lists of the
            entries from the supplied iterator.  Each list will have
            batch_size entries, although the final list may be shorter.
            
            src: https://biopython.org/wiki/Split_large_file
            """
            entry = True  # Make sure we loop once
            while entry:
                batch = []
                while len(batch) < batch_size:
                    try:
                        entry = next(iterator)
                    except StopIteration:
                        entry = None
                    if entry is None:
                        break
                    batch.append(entry)
                if batch:
                    yield batch
        # Open the large fasta file and use batch_iterator to split the file into batches of params.chunk_size sequences.
        os.makedirs(output.split_dir, exist_ok=True)
        record_iter = SeqIO.parse(open(input.fasta_file), "fasta")
        for i, batch in enumerate(batch_iterator(record_iter, params.chunk_size)):
            # Write the current batch to a split fasta file.
            output_file = f"{output.split_dir}/{i + 1}.fasta"
            handle = open(output_file, "w")
            SeqIO.write(batch, handle, "fasta")
            handle.close()
        with open(output.marker, "w") as f: # touches the file. Might come useful for storing information later
            f.write('')

import glob

def get_fasta_splits(wildcards):
    chkpt_output = checkpoints.split_fasta.get(**wildcards).output[0]
    return glob.glob(f"{chkpt_output}/*.fasta")

rule run_signalp:
    input: 
        fasta_file = "split_files/{file_id}.fasta",
        marker = "split_fasta.done"
    output:
        outfile = "split_sigp/{file_id}_summary.signalp5"
    params:
        prefix = "split_sigp/{file_id}"
    shell:
        """
        mkdir -p split_sigp
        signalp -batch 5000 -fasta {input.fasta_file} -org euk -format short -verbose -prefix {params.prefix}
        """

rule filter_signalp_outputs:
#   Description: this rule filters the output of the multiple signalp runs and extracts only those with a probability of signal peptide greater than a threshold. Only one file should be produced from the multiple signalp files. Two outputs are expected: a filtered (or not?) table with the signalp results and a filtered fasta of only those peptides with a signal
    input:
        files = expand("{file_id}", file_id = glob.glob("split_sigp/*.signalp5"))
    output:
        "filtered_sigp.tsv"
    params:
        threshold = 0.8 # todo: user defined
    shell:
        """
        cat {input.files} | sed '/^#/d' | awk '$3 > {params.threshold}' > {output}
        """

rule extract_secreted_peptides:
    input:
        signalp_result = rules.filter_signalp_outputs.output,
        fasta_file = rules.cluster_peptides.output.filtered_aa_sequences
    output:
        secreted_peptides = "secreted_peptides.fasta",
        non_secreted_peptides = "non_secreted_peptides.fasta"
    run:
        from Bio import SeqIO
        with open(str(input.signalp_result)) as infile:
            records = []
            for line in infile:
                records.append(line.rstrip().split("\t")[0])
        with open(output.non_secreted_peptides, "w") as n_outfile:
            with open(output.secreted_peptides, "w") as outfile:
                for seq in SeqIO.parse(input.fasta_file, "fasta"):
                    if seq.id in records:
                        SeqIO.write(seq, outfile, "fasta")
                    else:
                        SeqIO.write(seq, n_outfile, "fasta")



rule run_phobius: #todo: remember to inform the user about the installation procedure. I added a dependency in the conda env with a convenient installation script
#   Description: runs phobius
    input: 
        rules.extract_secreted_peptides.output.secreted_peptides
    output:
        table = "phobius_predictions.tsv"
    shell:
        """
        phobius -short {input} | sed 's/\s\+/\t/g' > {output.table}
        """

#todo: try to run signalp during the split rule to avoid problems. issue: if the process is interrupted abnormally during the run the rule is almost certain to misbehave and rerun the whole thing

#todo: rule run_signalp: # requires some testing. 
#todo: test with stripped sequences. This means that all sequences are preprocessed to be cut to a fixed length that would contain a signal peptide (like 50 or so). this might save memory and improve time. Moreover, we could try deduplicating these cut sequences and rereplicate afterwards to avoid predicting the same signal over and over. 
#   Description: runs signalp on the detected orfs

rule all: #todo: there is a bug that makes this rule run even if no split file is done. might solve by adding a checkpoint file at the end of the split rule
    input: 
        rules.extract_secreted_peptides.output,
        split_done = "split_fasta.done",
        signalp_files = expand("split_sigp/{f}_summary.signalp5", f=[i.split("/")[1].split(".")[0] for i in  glob.glob("split_files/*.fasta")])
