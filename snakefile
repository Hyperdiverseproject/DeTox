from Bio import SeqIO
import pandas 
import glob

def get_signalp_splits(wildcards):
    import os
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    return expand(global_output("")+'split_files/{i}.fasta',
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.fasta')).i)

def fastaToDataframe(fastaPath):
    sequences = SeqIO.parse(open(fastaPath), 'fasta')
    data = []
    for record in sequences:
        data.append({'ID': record.id, 'Sequence': str(record.seq)})
    return pandas.DataFrame(data)


def global_output(path):
    """ Construit le chemin d'output complet """
    output_dir = config.get("output_dir", "")
    if output_dir:
        return f"{output_dir}/{path}"
    else:
        return path

#configfile: "config.yaml"

if config['R1'] not in [None, ""]:
    if 'R2' in config and config['R2'] not in [None, ""]:
        # Reads paired-end
        rule trim_reads:
        #   Description: trims the provided raw reads
            input:
                r1 = config['R1'],
                r2 = config['R2'],
                adapters = config['adapters'],
            output:
                r1 = global_output("trimmed_reads/" + config['R1'].split("/")[-1]),
                r1_unpaired = global_output("trimmed_reads/unpaired." + config['R1'].split("/")[-1]),
                r2 = global_output("trimmed_reads/" + config['R2'].split("/")[-1]),
                r2_unpaired = global_output("trimmed_reads/unpaired." + config['R2'].split("/")[-1]),
            threads: config['threads']
            shell:
                """
                mkdir -p trimmed_reads
                trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} ILLUMINACLIP:{input.adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15
                """
    else:
        # Reads single-end
        rule trim_reads:
            # Description: trims the provided raw single-end reads
            input:
                r1 = config['R1'],
                adapters = config['adapters'],
            output:
                r1 = global_output("trimmed_reads/" + config['R1'].split("/")[-1]),
            threads: config['threads']
            shell:
                """
                mkdir -p trimmed_reads
                trimmomatic SE -threads {threads} {input.r1} {output.r1} ILLUMINACLIP:{input.adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15
                """

if config['transcriptome'] in [None, ""]:
    rule assemble_transcriptome:
        # Description: Assembles a transcriptome if it is not provided. Uses Trinity 
        # In this case sequencing reads MUST be provided in the config.
        input:
            r1 = rules.trim_reads.output.r1,
            r2 = lambda wildcards: rules.trim_reads.output.r2 if 'R2' in config and config['R2'] not in [None, ""] else [],
        output:
            assembly = global_output("trinity_out_dir/Trinity.fasta")
        params:
            memory = str(config['memory']) + "G",
            seqType = "fq",
            reads = lambda wildcards, input: "--single " + str(input.r1) if len(input.r2) == 0 else "--left " + str(input.r1) + " --right " + str(input.r2),
        threads: config['threads']
        shell:
            """
            Trinity --seqType {params.seqType} {params.reads} --CPU {threads} --max_memory {params.memory} --output {output.assembly}
            """

    
if 'contaminants' in config and config['contaminants'] not in [None, ""]:
    rule build_contaminants_database:
    #   Description: builds blast database for the removal of contaminants   
    #   todo: make this optional like in the original code
        input:
            fasta_db = config['contaminants']
        output:
            blast_db = global_output(config['contaminants'].split("/")[-1]+".out")
        shell:
            """
            touch {output.blast_db}
            makeblastdb -dbtype nucl -in {input.fasta_db} -out {output.blast_db}
            """

    rule blast_on_contaminants:
    #   Description: performs the actual blast of the contigs against the contaminants database
        input:
            blast_db = rules.build_contaminants_database.output.blast_db,
            contigs = lambda wildcards: config['transcriptome'] if 'transcriptome' in config and config['transcriptome'] not in [None, ""] else rules.assemble_transcriptome.output.assembly,
        output:
            blast_result = global_output(config['basename'] + ".blastsnuc.out")
        params:
            evalue = config['contamination_evalue']
        threads: config['threads']
        shell:
            """
            blastn -db {input.blast_db} -query {input.contigs} -out {output.blast_result} -outfmt 6 -evalue {params.evalue} -max_target_seqs 1 -num_threads {threads}
            """

    rule filter_contaminants:
    #   Description: performs the actual filtering
        input: 
            contigs = lambda wildcards: config['transcriptome'] if 'transcriptome' in config and config['transcriptome'] not in [None, ""] else rules.assemble_transcriptome.output.assembly,
            blast_result = rules.blast_on_contaminants.output.blast_result,    
        output:
            filtered_contigs = global_output(config['basename'] + ".filtered.fasta")
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
        nucleotide_sequences = rules.filter_contaminants.output.filtered_contigs if config['contaminants'] not in [None, ""] else config['transcriptome'] if 'transcriptome' in config  and config['transcriptome'] != None else rules.assemble_transcriptome.output.assembly
    output:
        aa_sequences = global_output(config['basename'] + ".faa")
    params:
        minlen = "99" if "minlen" not in config else config['minlen'],
        maxlen = "30000000" if "maxlen" not in config else config['maxlen']
    threads: config['threads']
    shell:
        """
        orfipy --procs {threads} --start ATG --partial-3 --partial-5 --pep {output.aa_sequences} --min {params.minlen} --max {params.maxlen} {input.nucleotide_sequences} --outdir .
        """

rule drop_X:
    input:
        aa_sequences = rules.detect_orfs.output.aa_sequences
    output:
        drop_sequence = global_output(config['basename'] + "_noX.faa")
    run:
        from Bio import SeqIO
        with open(f"{output}", "w") as outfile:
            for seq in SeqIO.parse(f"{input}", "fasta"):
                if "X" not in seq.seq:
                    SeqIO.write(seq, outfile, "fasta")

rule cluster_peptides:
#   Description: runs cd-hit on predicted peptide to remove excess redundancy
    input:
        aa_sequences = rules.drop_X.output.drop_sequence
    output:
        filtered_aa_sequences = global_output(config['basename'] + ".clustered.faa")
    params:
        threshold = config['clustering_threshold'], #todo : we might want to use separate thresholds if we're going to run cd-hit on transcripts and peptides
        memory = str(int(config['memory'])*1000),
        basename = config['basename'] + ".clustered"
    threads: config['threads']
    shell:
        """
        cd-hit -i {input.aa_sequences} -o {output.filtered_aa_sequences} -c {params.threshold} -M {params.memory} -T {threads} -d 40
        """

rule trim_peptides:
#   Description: this rule trims all peptides to only the first 50 aminoacids, as they are the only useful part for signalp. This step improves load time.
    input:
        aa_sequences = rules.cluster_peptides.output.filtered_aa_sequences
    output:
        trimmed_sequences = global_output(config['basename'] + ".trimmed.faa"),
    threads: 
        config['threads']
    run:
        from Bio import SeqIO
        import subprocess
        with open(output.trimmed_sequences, "w") as outfile:
            for seq in SeqIO.parse(input.aa_sequences, "fasta"):
                outfile.write(f">{seq.id}\n{seq.seq[:50]}\n")


checkpoint split_fasta:
    input:
        fasta_file = rules.trim_peptides.output.trimmed_sequences
    output:
        split_dir = directory(global_output("split_files"))
    params:
        chunk_size = 5000
    threads:
        config['threads']
    shell:
        """
        mkdir -p {output.split_dir}
        seqkit split2 -s {params.chunk_size} -O {output.split_dir} --by-size-prefix ""  -j {threads} {input.fasta_file}
        """


rule run_signalp:
    input: 
        fasta_file = global_output("")+"split_files/{i}.faa",
    output:
        outfile = global_output("split_files/{i}_summary.signalp5")
    params:
        prefix = global_output("split_files/{i}"),
        signalp_path = config['signalp_path']
    shell:
        """
        signalp -batch 5000 -fasta {input.fasta_file} -org euk -format short -verbose -prefix {params.prefix}
        """


def aggregate_splits(wildcards):
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    return expand(global_output("")+"split_files/{i}_summary.signalp5",
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.faa")).i)


rule filter_signalp_outputs:
#   Description: this rule filters the output of the multiple signalp runs and extracts only those with a probability of signal peptide greater than a threshold. Only one file should be produced from the multiple signalp files. Two outputs are expected: a filtered (or not?) table with the signalp results and a filtered fasta of only those peptides with a signal
    input:
        files = aggregate_splits
    output:
        global_output(config["basename"] + "_filtered_sigp.tsv")
    params:
        threshold = config['signalp_dvalue'] if 'signalp_dvalue' in config else "0.7"
    shell:
        """
        awk -v b={params.threshold} -F'\t' '!/^#/ && !/\?/  && $3 > b' {input.files} > {output}
        """

rule extract_secreted_peptides:
    input:
        signalp_result = rules.filter_signalp_outputs.output,
        fasta_file = rules.cluster_peptides.output.filtered_aa_sequences
    output:
        secreted_peptides = global_output(config['basename'] + "_secreted_peptides.fasta"),
        non_secreted_peptides = global_output(config['basename'] + "_non_secreted_peptides.fasta")
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
        table = global_output(config['basename'] + "_phobius_predictions.tsv")
    shell:
        """
        phobius.pl -short {input} | sed 's/\s\+/\t/g' | awk '$2 == 0' > {output.table}
        """

rule extract_non_TM_peptides:
#   Description: extracts non-TM peptides from the phobius output
    input:
        phobius_result = rules.run_phobius.output.table,
        fasta_file = rules.extract_secreted_peptides.output.secreted_peptides
    output:
        non_TM_peptides = global_output(config['basename'] + "non_TM_peptides.fasta")
    run:
        from Bio import SeqIO
        with open(str(input.phobius_result)) as infile:
            records = []
            for line in infile:
                records.append(line.rstrip().split("\t")[0])
        with open(output.non_TM_peptides, "w") as outfile:
            for seq in SeqIO.parse(input.fasta_file, "fasta"):
                if seq.id in records:
                    SeqIO.write(seq, outfile, "fasta")


rule build_toxin_blast_db: #todo: do we switch to diamond for max speed?
#   Description: builds a blast database for toxin prediction
    input:
        db = config['toxin_db']
    output:
        outfile = global_output(config['toxin_db'].split("/")[-1] + ".dmnd"),
    shell:
        """
        diamond makedb --db {output.outfile} --in {input.db} 
        """

rule blast_on_toxins:
#   Description: runs blastp against the toxin database. The query are the peptides without any signal sequence. The rule runs blast and extracts the fasta at the same time. might be split in two rules for easier management.
    input:
        orf_fasta_clustered_file = rules.cluster_peptides.output.filtered_aa_sequences,
        db_file = rules.build_toxin_blast_db.output.outfile,
    output:
        blast_result = global_output(config['basename'] + "_toxin_blast_results.tsv"),
    params:
        evalue = config['toxins_evalue'] if 'toxins_evalue' in config else "1e-10"
    threads: 
        config['threads']
    run:
        import subprocess
        from Bio import SeqIO
        build_header = f"echo \"qseqid\ttoxinDB_sseqid\ttoxinDB_pident\ttoxinDB_evalue\" > {output.blast_result}"
        command_line = f"{build_header} && diamond blastp -q {input.orf_fasta_clustered_file} --evalue {params.evalue} --max-target-seqs 1 --threads {threads} -d {input.db_file} --outfmt 6 qseqid sseqid pident evalue >> {output.blast_result}"
        print(command_line)
        subprocess.run(command_line, shell=True)
        

rule retrieve_orfs_with_blast_without_signalp:
    input:
        nonsec_fasta_file = rules.extract_secreted_peptides.output.non_secreted_peptides,
        blast_toxins_result = rules.blast_on_toxins.output.blast_result
    output:
        hits_fasta = global_output(config['basename'] + '_toxins_by_similarity.fasta')
    run:
        with open(str(input.blast_toxins_result)) as infile:
            records = []
            for line in infile:
                records.append(line.rstrip().split("\t")[0])
        with open(output.hits_fasta, "w") as outfile:
            for seq in SeqIO.parse(input.nonsec_fasta_file, "fasta"):
                if seq.id in records:
                    SeqIO.write(seq, outfile, "fasta")


rule retrieve_candidate_toxins:
#   Description: this rule just creates a fasta from the positive hits in the toxin similarity and structure searches. 
    input:
        structure_based = rules.extract_non_TM_peptides.output.non_TM_peptides,
        similarity_based = rules.retrieve_orfs_with_blast_without_signalp.output.hits_fasta
    output:
        global_output(config["basename"] + "_candidate_toxins.fasta")
    shell:
        """
        cat {input.structure_based} {input.similarity_based} > {output}
        """


rule download_pfam:
    output:
        db_dir = directory(global_output("databases/pfam")),
        pfam_db = global_output("databases/pfam/Pfam-A.hmm")
    shell:
        """
        mkdir -p {output.db_dir}
        wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -P {output.db_dir}
        gunzip -c {output.db_dir}/Pfam-A.hmm.gz > {output.pfam_db}
        """

rule run_hmmer:
#   Description: runs hmmer against the pfam database. 
    input:
        fasta_file = rules.retrieve_candidate_toxins.output,
        pfam_db = lambda wildcards: config["pfam_db_path"] if "pfam_db_path" in config and config["pfam_db_path"] not in [None, ""] else rules.download_pfam.output.pfam_db,
    output:
        tblout = global_output(config['basename'] + ".tblout"),
        domtblout = global_output(config['basename'] + ".domtblout")
    params:
        evalue = config['pfam_evalue'] if 'pfam_evalue' in config else "1e-5"
    threads: 
        config['threads']
    shell:
        """
        hmmsearch --cut_ga --cpu {threads} --domtblout {output.domtblout} --tblout {output.tblout} {input.pfam_db} {input.fasta_file} 
        """

rule parse_hmmsearch_output:
#   Description: parses and aggregates the hmmer output, uses the domtblout file
    input: 
        domtblout = rules.run_hmmer.output.domtblout
    output:
        filtered_table = global_output(config['basename'] + ".domtblout.tsv")
    run:
        import pandas
        df_domtblout = pandas.read_csv(f"{input.domtblout}", comment="#", delim_whitespace=True, usecols = [0,1,2,3,4] , names=["target name","accession_t","tlen","query name","accession_Q"])
        aggregated_domains = df_domtblout.groupby('target name')['query name'].apply(list).reset_index()
        aggregated_domains['pfam domains'] = aggregated_domains['query name'].apply(lambda x: "; ".join(list(set(x))))
        aggregated_domains.to_csv(f"{output.filtered_table}", sep="\t", index=False)

rule run_wolfpsort:
#   Description: runs wolfpsort on secreted peptides inferred by signalp 
    input:
        rules.retrieve_candidate_toxins.output
    output:
        global_output(config['basename'] + "_secreted_wolfpsort_prediction.tsv")
    params:
        wps_path = config['wolfPsort_path'],
        awk = "awk '{print $1\"\t\"$2}'"
    shell:
        """
        {params.wps_path} animal < {input} | {params.awk} > {output}
        """

rule detect_repeated_aa:
#   Description: this rule looks at the fasta aminoacid sequences in input and produces a table. The table reports whether some kind of repeated pattern is found in the sequences (up to 3AA long). The default threshold for repetition is 5. The input is processed with biopython
    input:
        fasta_file = rules.retrieve_candidate_toxins.output,
    output:
        repeated_aa = global_output(config['basename'] + "_repeated_aa.tsv")
    params:
        threshold = config['repeated_aa_threshold'] if 'repeated_aa_threshold' in config else 5
    run:
        import itertools
        from Bio import SeqIO
        import pandas
        def findRepetition(size,seq):
            repetition=[]
            for cdl in range(0,size):
                sub = [seq[i:i+size] for i in range(cdl, len(seq), size)]
                groups = itertools.groupby(sub)
                result = [(label, sum(1 for _ in group)) for label, group in groups]
                for elem, nbRep in result:
                    if int(nbRep) >=int(f"{params.threshold}"):
                        repetition.append((elem,nbRep))
            return repetition
        secreted = fastaToDataframe(f"{input.fasta_file}")
        secreted["Repeats1"] = secreted.apply(lambda x: findRepetition(1,x["Sequence"]),axis=1)
        secreted["Repeats2"] = secreted.apply(lambda x: findRepetition(2,x["Sequence"]),axis=1)
        secreted["Repeats3"] = secreted.apply(lambda x: findRepetition(3,x["Sequence"]),axis=1)
        secreted["Repeats"] = secreted["Repeats1"] + secreted["Repeats2"] + secreted["Repeats3"]
        secreted['RepeatsTypes'] = secreted['Repeats'].apply(lambda t: [n for (n, _) in t])
        secreted['RepeatsLengths'] = secreted['Repeats'].apply(lambda t: [n for (_, n) in t])
        secreted['RepeatsLengths'] = [','.join(map(str, l)) for l in secreted['RepeatsLengths']]
        secreted['RepeatsTypes'] = [','.join(map(str, l)) for l in secreted['RepeatsTypes']]
        secreted = secreted.drop(columns=["Repeats","Repeats1","Repeats2","Repeats3"])
        secreted.to_csv(f"{output.repeated_aa}", index=False,sep='\t')

def get_cys_pattern(seq):
    pattern = ""
    status = False
    if not pandas.isna(seq) and seq.count('C') >= 4:
        for char in seq:
            if char == "C":
                pattern = pattern + "C"
                status = True
            else:
                if status:
                    pattern = pattern + "-"
                    status = False
        if(pattern[-1] == "-"):
            pattern = pattern[0:-1]
    if pattern == "":
        pattern = None
    return pattern


### optional rules

if config['quant'] == True:
    if config['R1'] not in [None, ""]:
        if 'R2' in config and config['R2'] not in [None, ""]:
            # Reads paired-end
            rule run_salmon:
                input:
                    transcriptome = lambda wildcards: config['transcriptome'] if 'transcriptome' in config and config['transcriptome'] not in [None, ""] else rules.assemble_transcriptome.output.assembly,
                    r1 = rules.trim_reads.output.r1,
                    r2 = rules.trim_reads.output.r2
                output:
                    quant_dir = directory(global_output(config['basename'] + "_quant")),
                    index= directory(global_output("salmon.idx")),
                    quantification = global_output(config['basename'] + "_quant/quant.sf")
                threads:
                    config['threads']
                shell:
                    """
                    salmon index -t {input.transcriptome} -i {output.index} -p {threads}
                    salmon quant -i {output.index} -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {output.quant_dir}
                    """
        else:
            rule run_salmon:
                input:
                    transcriptome = lambda wildcards: config['transcriptome'] if 'transcriptome' in config and config['transcriptome'] not in [None, ""] else rules.assemble_transcriptome.output.assembly,
                    r1 = rules.trim_reads.output.r1
                output:
                    quant_dir = directory(global_output(config['basename'] + "_quant")),
                    index= directory(global_output("salmon.idx")),
                    quantification = global_output(config['basename'] + "_quant/quant.sf")
                threads:
                    config['threads']
                shell:
                    """
                    salmon index -t {input.transcriptome} -i {output.index} -p {threads}
                    salmon quant -i {output.index} -l A -r {input.r1} --validateMappings -o {output.quant_dir}
                    """


rule download_uniprot:
    output:
        db_dir = directory(global_output("databases/uniprot")),
        database = global_output("databases/uniprot/uniprot_sprot.fasta.gz")
    shell:
        """
        mkdir -p {output.db_dir}
        wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -P {output.db_dir}
        """


rule make_uniprot_blast_database:
#   Description: builds a blast database from the uniprot fasta
    input:
        fasta_file = lambda wildcards: config["swissprot_db_path"] if "swissprot_db_path" in config and config["swissprot_db_path"] not in [None, ""] else rules.download_uniprot.output.database
    output:
        db_file = global_output("uniprot_blast_db.dmnd")
    shell:
        """
        diamond makedb --in {input.fasta_file} --db {output.db_file}
        """

rule blast_on_uniprot:
#   Description: run blast against the uniprot database, return only the best hit
    input:
        fasta_file = rules.retrieve_candidate_toxins.output,
        db_file = rules.make_uniprot_blast_database.output.db_file
    output:
        blast_result = global_output(config['basename'] + "_uniprot_blast_results.tsv")
    params:
        evalue = config['swissprot_evalue'] if 'swissprot_evalue' in config else "1e-10",
    threads: 
        config['threads']
    shell:
        """
        echo "qseqid\tuniprot_sseqid\tuniprot_pident\tuniprot_evalue" > {output.blast_result}
        diamond blastp -d {input.db_file} -q {input.fasta_file} --evalue {params.evalue} --outfmt 6 qseqid sseqid pident evalue --max-target-seqs 1 --threads {threads} >> {output.blast_result}
        """


# TODO: follow this comment for the rule that will wraps everything up and create the final table. -> Also, in my opinion these peptides should be marked with a warning flag in the output, specifying which issue affects them (e.g. “this peptide lacks a signal peptide”, “this peptide contains a transmembrane domain”, etc.)



#todo: try to run signalp during the split rule to avoid problems. issue: if the process is interrupted abnormally during the run the rule is almost certain to misbehave and rerun the whole thing

# this is the list with all the expected output to be put in the final table, will be filled depending on the configuration file. 
outputs = [
    (rules.run_wolfpsort.output if config['wolfpsort'] else []),
    rules.parse_hmmsearch_output.output,
    rules.blast_on_toxins.output.blast_result,
    (rules.blast_on_uniprot.output.blast_result if config['swissprot'] else []),
    rules.detect_repeated_aa.output.repeated_aa,
]


rule build_output_table:
#   Description: this rule merges the tabular output of the other rules and merges it in a single table. It uses the outputs list defined above.
    input:
        #base = rules.extract_Cys_pattern.output,
        fasta_file = rules.retrieve_candidate_toxins.output,
        signalp_result = rules.filter_signalp_outputs.output,
        quant = rules.run_salmon.output.quantification if config['quant'] else [],
        extra = outputs
    output:
        global_output(config['basename'] + "_toxins.tsv")
    params:
        TPMthreshold = config['TPMthreshold'] if 'TPMthreshold' in config else 1000,
    run:
        seq_df = fastaToDataframe(f"{input.fasta_file}")
        signalp_df = pandas.read_csv(f"{input.signalp_result}", sep='\t', names = ["ID", "signalp_prediction", "prob_signal", "prob_nosignal", "cutsite"])
        df = seq_df.merge(signalp_df, on = "ID", how = "left")
        #df = pandas.read_csv(f"{input.base}", sep='\t')
        for i in input.extra:
            dfi = pandas.read_csv(str(i), sep = "\t")
            if "Sequence" in dfi.columns:
                dfi = dfi.drop(columns=["Sequence"])
            newcols = [i for i in dfi.columns]
            newcols[0] = "ID"
            dfi.columns = newcols
            df = df.merge(dfi, how = "left", on = "ID")
        if config['cys_pattern'] == True:
            df['cutsite'] = df['cutsite'].fillna("")
            df['cut_site_position'] = df['cutsite'].apply(lambda x: int(x.split(" ")[2].split("-")[-1][:-1]) if "pos:" in x else -1)
            df['mature_peptide'] = df.apply(lambda x: x['Sequence'][x['cut_site_position']:] if x['cut_site_position']>0 else None , axis=1)
            df['Cys_pattern'] = df['mature_peptide'].apply(lambda x: get_cys_pattern(x) if pandas.notna(x) else None)
        if config['quant'] == True:
            df['contig'] = df['ID'].apply(lambda x: x.split("_ORF")[0])
            q = pandas.read_csv(f"{input.quant}", sep = "\t")
            newcols = [i for i in q.columns]
            newcols[0] = "contig"
            q.columns = newcols
            df = df.merge(q, how = "left", on = "contig")
            df = df.drop(['EffectiveLength','NumReads'], axis=1)
        try:
            df = df.assign(Rating='')
            df['Rating'] = df.apply(lambda row: str(row['Rating'] + 'S') if pandas.notna(row['signalp_prediction']) else str(row['Rating'] + '*'), axis=1)
            df['Rating'] = df.apply(lambda row: str(row['Rating'] + 'B') if pandas.notna(row['toxinDB_sseqid']) else row['Rating'], axis=1)
            if 'Cys_pattern' in df.columns:
                print("cys pattern et son type : ",df["Cys_pattern"][0],type(df["Cys_pattern"][0]))
                df['Rating'] = df.apply(lambda row: str(row['Rating'] + 'C') if pandas.notna(row['Cys_pattern']) else row['Rating'], axis=1)
            if 'TPM' in df.columns:
                df['Rating'] = df.apply(lambda row: str(row['Rating'] + 'T') if (float(row['TPM'])>=float(f"{params.TPMthreshold}")) else row['Rating'], axis=1)
            df['Rating'] = df.apply(lambda row: str(row['Rating'] + 'D') if pandas.notna(row['pfam domains']) else row['Rating'], axis=1)
            if 'uniprot_sseqid' in df.columns:
                df['Rating'] = df.apply(lambda row: str(row['Rating'] + '!') if pandas.notna(row['uniprot_sseqid']) and pandas.isna(row['toxinDB_sseqid']) else row['Rating'], axis=1)
        except Exception as e:
            print(f"An error has occurred during sequence rating: {e}")
            sys.exit()  
        #df = df[df['mature_peptide'].apply(lambda x: len(str(x))) > 3]
        df = df.drop(['cut_site_position','query name'], axis=1)
        df.rename(columns={'k': 'wolfpsort_prediction'}, inplace=True)
        df.drop_duplicates().to_csv(f"{output}", sep='\t', index=False)

rule all:
    input: 
       rules.build_output_table.output
