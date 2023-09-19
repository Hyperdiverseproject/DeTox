from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import subprocess
import multiprocessing as mp
import os
from Bio import SeqIO
import urllib
import gzip
import shutil
import textwrap
import pandas
from io import StringIO
import sys
import itertools
from shutil import which
from Bio.Blast.Applications import NcbiblastpCommandline

#================================ RECUPERATION DES PARAMETRES ======================================================================================================================================

print("|-----------------------------------------------------|")
print("|                    DeTox                            |")
print("|-----------------------------------------------------|")

parser = ArgumentParser(prog='PROG',
         formatter_class=RawDescriptionHelpFormatter,
         description=textwrap.dedent('''\
             All Programs must be installed !
            --------------------------------
            It is preferable to install the programs with "sudo".
            --------------------------------
            
            Here is the list of programs used :
                - ORFfinder
                - Signalp
                - Phobius
                - hmmsearch
                - makeblastdb
                - blastp
                - cd-hit
                
                OPTIONAL :
                
                    With assembly:
                
                        - Trinity  needs : Jellyfish, Salmon, samtools, Bowtie
                        - Trimmomatic 
                        
                    Other filters:
                    
                        - Salmon (opt)
                        - RSEM  (opt) ( long to run )
                        - WoLF PSORT  (opt)
            
            --------------------------------
            
            Here is the list of files needed for the program :
                - Transcriptome (fasta)
                - Reads files of the transcriptome ( fastq )
                - Toxin Database (fasta)
                
             '''))

parser.add_argument("-A",
                      help="Create the assembly with Trinity and its controls from the reads in fastq",
                      dest="Assembly",
                      action='store_true')

parser.add_argument("-M",
                      type=str,
                      help="gigabit available memory",
                      dest="memory",
                      default="8")

parser.add_argument("-t",
                      type=str,
                      help="Path of Transcriptome, file name format: specie.fasta",
                      dest="TranscFilesPath")

parser.add_argument("-S",
                      help="Adds Salmon quantification",
                      dest="Salmon",
                      action='store_true')

parser.add_argument("-r1",
                      type=str,
                      help="Path of R1 reads file",
                      dest="R1ReadsFile")

parser.add_argument("-r2",
                      type=str,
                      help="Path of R2 reads file",
                      dest="R2ReadsFile")

parser.add_argument("-dbT",
                      type=str,
                      help="Path of toxin database",
                      dest="ToxinDatabase",
                      required=True)

parser.add_argument("-dbC",
                      type=str,
                      help="Path of contamination database",
                      dest="ContaminationDatabase")

parser.add_argument("-dbP",
                      help="Blast on Swiss-Prot database, it will be downloaded if no path is provided",
                      dest="SwissProtDatabase",
                      action='store_true')

parser.add_argument("-W",
                      type=str,
                      help="Path to WolpSort program 'runWolfPsortSummary' ",
                      dest="wolfPSortPath")

parser.add_argument("-cys",
                      help="Run cys pattern detection",
                      dest="cysPattern",
                      action='store_true')

parser.add_argument("-EC",
                      type=str,
                      help="Contamination Evalue threshold",
                      dest="EvalueContamination",
                      default="1E-5")

parser.add_argument("-C",
                      type=str,
                      help="CD-Hit clustering threshold",
                      dest="ClusteringThreshold",
                      default="0.99")

parser.add_argument("-Ol",
                      type=str,
                      help="ORFfinder minimum length of the reading frame",
                      dest="ORFlength",
                      default="33")

parser.add_argument("-tSP",
                      type=str,
                      help=" only sequences with a signal peptide having a probability greater than the threshold are kept",
                      dest="SignalPeptideProbability",
                      default="0.7")

parser.add_argument("-ET",
                      type=str,
                      help="Blastp Evalue on the toxin database",
                      dest="BlastpToxinesEvalue",
                      default="1E-10")

parser.add_argument("-R",
                      type=str,
                      help="Minimal length of amino acid repeat to be consider",
                      dest="RepeatsLength",
                      default="5")

parser.add_argument("-ES",
                      type=str,
                      help="Blastp Evalue on the SwissProt database",
                      dest="BlastpSwissprotEvalue",
                      default="1E-5")
parser.add_argument("-tpm",
                      type=str,
                      help="TPM threshold used to select highly expressed sequences",
                      dest="TPMthreshold",
                      default="1000")

args = parser.parse_args()

#================================ Check parameters ===================================================================================================================================================

if (args.Salmon) and (args.R1ReadsFile is None or args.R2ReadsFile is None):
    parser.error("-S requires -r1 and -r2")
    
# The assembly needs to specify the memory
if args.Assembly and args.memory is None:
    parser.error("-A requires -M")
    
# We have to give a transcriptome or make an assembly
if not args.Assembly and args.TranscFilesPath is None:
    parser.error("-A or -t require")

name = ""
transcFilesPath = ""
# Get the name of the specie from the transcriptome path
if args.TranscFilesPath is not None:
    transcFilesPath = args.TranscFilesPath
    name = args.TranscFilesPath.split("/")[-1].split(".")[0]

finalStat = dict()

#================================ Check Programs ===================================================================================================================================================

listProgram = ["cd-hit","trimmomatic","Trinity","makeblastdb","blastn","ORFfinder","signalp","phobius","hmmsearch"]
FNULL = open(os.devnull, 'w')
programError = []

for program in listProgram:
    if which(program) is not None:
        print ("{:<20} | --------------------------------------------------| OK".format(program))
    else:
        print ("{:<20} | --------------------------------------------------| ERROR".format(program))
        programError.append(program)

if len(programError)>0:
    print("Les programmes suivant ne sont pas installer ou ne sont pas accessible: ",[elem for elem in programError])
    sys.exit()

#=============================== Function =======================================================================================================================================================================

def fastaToDataframe(fastaPath):
    sequences = SeqIO.parse(open(fastaPath), 'fasta')
    data = []
    for record in sequences:
        data.append({'ID': record.id, 'Sequence': str(record.seq)})
    return pandas.DataFrame(data)

def dataframe_to_fasta(df, sequence_column_name, id_column_name,outputFile_name):
    fasta = ""
    for index, row in df.iterrows():
        fasta += ">{}\n{}\n".format(row[id_column_name], row[sequence_column_name])
    with open(outputFile_name, 'w') as f:
        f.write(fasta)


#================================ Filtering READS QUALITY BEFORE ASSEMBLY ===================================================================================================================================================

if args.Assembly :
    nameR1reads = args.R1ReadsFile.split("/")[-1].split(".")[0]
    nameR2reads = args.R2ReadsFile.split("/")[-1].split(".")[0]
    try:
        subprocess.Popen(["trimmomatic PE "+args.R1ReadsFile+" "+args.R2ReadsFile+" "+nameR1reads+"_paired.fq.gz "+nameR1reads+"_unpaired.fq.gz "+nameR2reads+"_paired.fq.gz "+nameR2reads+"_unpaired.fq.gz ILLUMINACLIP:/home/allan/Program/Trimmomatic-0.39/adapters/adapters.fa:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15"],shell=True).communicate()
    except FileNotFoundError:
        print("TrimmomaticPE is not installed or was not found by the pipeline")
        sys.exit() 

#================================ ASSEMBLY WITH TRINITY ===================================================================================================================================================
if args.Assembly :
    CPU = str(mp.cpu_count())
    try:
        subprocess.Popen(["Trinity","--seqType","fq","--left",(nameR1reads+"_paired.fq.gz"),"--right",(nameR2reads+"_paired.fq.gz"),"--CPU",CPU,"--max_memory",(args.memory+"G"),"--KMER_SIZE","31"]).communicate()
        transcFilesPath = "trinity_out_dir/Trinity.fasta"
        name = transcFilesPath.split("/")[-1].split(".")[0]
    except FileNotFoundError:
        print("Trinity is not installed or was not found by the pipeline")
        sys.exit()



#================================ Filtering contaminations from the assembly ===================================================================================================================================================

if args.ContaminationDatabase is not None:
    
    def readBlastConta(filename):
        records = []
        infile = open(filename, 'r')
        for line in infile:
            line = line.rstrip()
            if line[0] != '#':
                blast = line.split()                                        
                records.append(blast[0]) # On recupère l'ID
        infile.close()
        return records
    
    def filterContaInTranscrit(listIDConta,transcitFile):
        recordIter = SeqIO.parse(open(transcitFile), "fasta")
        nbNonConta = 0
        nbConta =0
        with open(name+"_filtered.fasta", "w") as handle:
            for rec in recordIter:
                if rec.id not in listIDConta:
                    SeqIO.write(rec, handle, "fasta")
                    nbNonConta+=1
                else:
                    nbConta+=1
        return (nbNonConta,nbConta)
    
    try:
        subprocess.Popen(["makeblastdb","-in",args.ContaminationDatabase,"-dbtype","nucl"]).communicate()
        subprocess.Popen(["blastn","-db",args.ContaminationDatabase,"-query",transcFilesPath,"-out",name+".blastsnuc.out","-outfmt","6","-evalue",args.EvalueContamination,"-max_target_seqs","1","-num_threads",str(mp.cpu_count())]).communicate()        
        sequenceConta=readBlastConta(name+".blastsnuc.out")
        statsConta = filterContaInTranscrit(sequenceConta,transcFilesPath)
        finalStat['NB Seq with contamination :'] = statsConta[1]
        finalStat['NB Seq without contamination :'] = statsConta[0]
        transcFilesPath=name+"_filtered.fasta"
    except FileNotFoundError:
        print("The blast program is not installed or was not found by the pipeline")
        sys.exit()

     
#================================ Clustering contigs of ASSEMBLY ===================================================================================================================================================

try:
    CPU = str(mp.cpu_count())
    memory = args.memory+"000"
    subprocess.Popen(["cd-hit","-i",transcFilesPath,"-o",(name+"Clustered.fasta"),"-c",args.ClusteringThreshold,"-M",memory,"-T",CPU]).communicate()
    transcFilesPath = (name+"Clustered.fasta")
except FileNotFoundError:
    print("cd-hit is not installed or was not found by the pipeline")
    sys.exit()

#================================ ORFs DETECTION ===================================================================================================================================================

try:
    subprocess.Popen(["ORFfinder","-in",transcFilesPath,"-s","0","-ml",args.ORFlength,"-out","CDS-detected.fa"]).communicate()
except FileNotFoundError:
    print("The ORFfinder program is not installed or was not found by the pipeline")
    sys.exit()


#================================ Clustering proteins of ORFS detection ===================================================================================================================================================

try:
    CPU = str(mp.cpu_count())
    memory = args.memory+"000"
    subprocess.Popen(["cd-hit","-i","CDS-detected.fa","-o","CDS-Clustered.fasta","-c","1","-M",memory,"-T",CPU]).communicate()
    transcFilesPath = (name+"Clustered.fasta")
except FileNotFoundError:
    print("cd-hit is not installed or was not found by the pipeline")
    sys.exit()

#================================ PREPARATION AND RUN OF SIGNALP ==================================================================================================================================

# Spécifiez le chemin d'accès à votre fichier fasta volumineux.
fasta_file = "./CDS-Clustered.fasta"
dirTest = f"./split_fasta_{name}"
seqs_per_file = 50000

if not os.path.exists(dirTest):
    os.makedirs(dirTest)

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

# Ouvrir le fichier fasta volumineux et utiliser batch_iterator pour fragmenter le fichier en lots de 50 000 séquences.
record_iter = SeqIO.parse(open(fasta_file), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, seqs_per_file)):
    # Écrire le lot actuel dans un fichier fasta fragmenté.
    output_file = f"{dirTest}/{i + 1}.fasta"
    SeqIO.write(batch, output_file, "fasta")

# Lister tous les fichiers fasta fragmentés dans le dossier "split_fasta".
files = os.listdir(dirTest)
fasta_files = [f for f in files if f.endswith(".fasta")]

# Définir la fonction à exécuter en parallèle.
def run_signalp(fasta_file):
    input_file = f"{dirTest}/{fasta_file}"
    output_file = f"{dirTest}/{fasta_file}"
    cmd = f"signalp -fasta {input_file} -org euk -format short -verbose -batch {seqs_per_file} -prefix {output_file}"
    subprocess.call(cmd.split())

# Exécutez la fonction run_signalp en parallèle sur chaque fichier fasta fragmenté.
with mp.Pool(processes=os.cpu_count()) as pool:
    pool.map(run_signalp, fasta_files)

#================================ MERGE ALL FILES ===================================================================================================================================================

subprocess.Popen([f"cat {dirTest}/*_summary.signalp5 > {name}.signalp"],shell=True).communicate()

subprocess.Popen([f"rm -r {dirTest}"],shell=True).communicate()

#================================ CHECKS RESULTS OF THE SignalP RUN AND COPY THEN IN NEW FILES =======================================================================================================

def checkSignalpResultAndParse(signalPFilePath,ORFfilePath,Dvalue,outORFsNameFile,outSSNameFile):
    """CHECKS RESULTS OF THE SignalP RUN, the file produce by signalp doesn't contain the sequence but only ID
    So this function use IDs of signalp output to get sequence in original ORFs file
    
    Le fichier de signalp contient toute les potentiel ORF détecter mais beaucoup sont classer comme "OTHER". Celle ci ne
    nous interesse pas car elles n'ont pas été identifier. Ce filtre permet de conserver uniquement les séquence identifier et qui
    on une DValue > a un seuil.

    Args:
        signalPFilePath (str): _description_
        ORFfilePath (str): _description_
        Dvalue (float): Probality value
        outORFsNameFile (str): sequence of ORFs containing a signal peptide
        outSSNameFile (str): signal sequence after cleavage
    """
    sigdata = open(signalPFilePath, 'r')  # On ouvre le fichier de sequences signalp
    signalin = []
    sigends = []
    for line in sigdata:
        line = line.strip()     # On supprime les caractères de début et de fin
        if line[0] !='#':       # si le premiere caractere de la ligne est différent de #      
            data = line.split('\t')         # On split la ligne sur les tabulations
            if data[1] != 'OTHER' and float(data[2]) >= float(Dvalue):         # si le deuxième éléments est différent de OTHER et que la probabilité d'avoir un signal peptide ( colonne 3 ) est supérieur a la DValue
                signalin.append(data[0])            # On ajoute l'identifiant de l'ORF        
                sigend = line[line.index('CS pos: ')+len('CS pos: '): line.index('CS pos: ')+len('CS pos: ')+2].split("-")[0] # coordonné du site de clivage dans la séquence
                if '?' in sigend:
                    sigends.append('q')             # ajouter 'q' si il y a un ? dans les coordonnées
                else:
                    sigends.append(int(sigend))       # ajouter les coordonnées du site de clivage dans sigends
    sigdata.close()                                     # On ferme le fichier 
    print('ORFs with signal sequence:', len(signalin), len(sigends))
    finalStat['ORFs with signal sequence:'] = len(signalin)
    print('ORFs with unresolved Signal Séquence boundary:', sigends.count('q'))
    finalStat['ORFs with unresolved Signal Séquence boundary:'] = sigends.count('q')
    nbORF=0
    allORFsfasta = open(ORFfilePath, 'r')        # on ouvre le fichier d'ORF
    outORFsfasta = open(outORFsNameFile, 'w')       # on ouvre le fichier d'écriture des séquence des ORFs contenant un peptide signal
    outSSfasta = open(outSSNameFile, 'w')           # fichier ou on écrit les signals sequences après clivage
    for line in allORFsfasta:               # pour chaque ligne des ORFs 
        line=line.strip()                       # on supprime les caractères d'espacement de début et de fin
        if line.startswith('>'):                # si la ligne commence par un '>'
            nbORF+=1
            ORFID = line[1:].split()[0]             # on recupere l'id de l'ORF qui se trouve en première position après avoir split sur l'espace
            if ORFID in signalin:                       # si cette ORF est presente dans la liste des ORF detecter par signalp
                seq = next(allORFsfasta).strip()       # On recupere la sequence de l'ORF
                ind = signalin.index(ORFID)             # on recupère l'indice de la sequence dans la liste des ORF avec peptide signal
                send = sigends[ind]                     # on recupére les coordonné de clivage dans la séquence ORF
                if send != 'q':                                 # si c'est différent de q 
                    outORFsfasta.write('>'+ORFID+'\n')                  # On réécrit la séquence complete
                    outORFsfasta.write(seq.replace('*', '')+'\n')       # de l'ORF dont on a détecter un peptide signal
                    
                    outSSfasta.write('>'+ORFID+'\n')                    # On écrit uniquement la séquence de la protéine 
                    outSSfasta.write(seq[:send]+'\n')                   # On recupere la séquence signal contenu entre le début et le site de clivage
    finalStat['NB ORFs detected by ORFfinder:'] = nbORF
    allORFsfasta.close()
    outORFsfasta.close()
    outSSfasta.close()

checkSignalpResultAndParse((name+".signalp"),"CDS-Clustered.fasta",args.SignalPeptideProbability,(name+"_sigORFs.fasta"),(name+"_SSeqs.fasta"))



#================================ DETECT TRANSMEMBRANE DOMAINS ======================================================================================================================================
print('phobius')

try:
    subprocess.Popen(['phobius -short '+name+"_sigORFs.fasta > "+name+'.phobius.out'],shell=True).communicate()
except FileNotFoundError:
    print("The phobius.pl program is not installed or was not found by the pipeline")
    sys.exit()

#================================ FILTERED ORF SEQUENCE OF SIGNALP WITH PHOBIUS OUTPUT ===============================================================================================================

phobius_df = pandas.read_csv((name+'.phobius.out'), sep='\s+', skiprows=1, names=['SequenceID', 'TM', 'SP', 'Prediction'])
signalpOutput_df= fastaToDataframe((name+'_sigORFs.fasta'))
phobius_filter_df = pandas.merge(phobius_df,signalpOutput_df,left_on="SequenceID",right_on="ID")
phobius_filter_df = phobius_filter_df[phobius_filter_df["TM"]==0]
phobius_filter_df = phobius_filter_df.drop(["TM","SP","Prediction","SequenceID"],axis=1, errors='ignore')
dataframe_to_fasta(phobius_filter_df, "Sequence", "ID",(name+'_finalORFs.fasta'))

#================================ Similarity approach, detection of toxins from ORFs with a blastp on a toxin database ================================================================================================================================================

try:
    subprocess.Popen(["makeblastdb","-in",args.ToxinDatabase,"-dbtype","prot"]).communicate()
    subprocess.Popen(["blastp","-db",args.ToxinDatabase,"-query",name+'_finalORFs.fasta',"-out",name+"_finalORFs.blastsprot.out","-outfmt","6","-evalue",args.BlastpToxinesEvalue,"-max_target_seqs","1","-num_threads","8"]).communicate()
except FileNotFoundError:
    print("The blast program is not installed or was not found by the pipeline")
    sys.exit()

#================================ Search of toxines in ORF without Signal peptide =====================================================================================================================================================================================

fasta_ORFs = "CDS-Clustered.fasta"
fasta_ORFs_withSP = name+'_finalORFs.fasta'
fasta_ORFs_withoutSP = name+"_SequenceWithoutSP.fa"
ORFS_withoutSP_with_blastHit = name+"_ORFS_withoutSP_with_blastHit.out"

# Liste des identifiants dans le fichier FASTA 2
ids_fasta2 = [record.id for record in SeqIO.parse(fasta_ORFs_withSP, "fasta")]

# Liste des séquences dans le fichier FASTA 1 mais pas dans le fichier FASTA 2
sequences_only_in_fasta1 = []
for record in SeqIO.parse(fasta_ORFs, "fasta"):
    if record.id.split(" ")[0] not in ids_fasta2:
        record.id=record.id.split(" ")[0]
        record.description=record.description.split(" ")[0]
        sequences_only_in_fasta1.append(record)

# Écrit les séquences dans le fichier de sortie
with open(fasta_ORFs_withoutSP, "w") as f:
    SeqIO.write(sequences_only_in_fasta1, f, "fasta")

dfSequenceWithoutSP = fastaToDataframe(fasta_ORFs_withoutSP)
blastp_cline = NcbiblastpCommandline(query=fasta_ORFs_withoutSP, db=args.ToxinDatabase, evalue=args.BlastpToxinesEvalue, outfmt=6,max_target_seqs=1,out=ORFS_withoutSP_with_blastHit)
blastp_cline()
blast_results_withoutSP = pandas.read_csv(ORFS_withoutSP_with_blastHit, sep='\t', header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
df_ORFs_withoutSP_whithToxinBLast = pandas.merge(dfSequenceWithoutSP,blast_results_withoutSP,left_on="ID",right_on="qseqid")

# ADD blast result of ORF without SP in blast file
blast_df_withSP = pandas.read_csv(name+"_finalORFs.blastsprot.out", sep="\t",names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
blast_df = pandas.concat([blast_df_withSP, blast_results_withoutSP], ignore_index=True)
fichier_blast = open(name+"_finalORFs.blastsprot.out", "w")
blast_df.to_csv(fichier_blast, sep="\t", index=False, header=False)
fichier_blast.close()

# ADD sequence of ORF without SP in finalORF file

dfSequenceWithSP = fastaToDataframe(fasta_ORFs_withSP)
df_ORFs_withoutSP_whithToxinBLast = df_ORFs_withoutSP_whithToxinBLast[["ID","Sequence"]]
df_All_final_ORF = pandas.concat([df_ORFs_withoutSP_whithToxinBLast, dfSequenceWithSP], ignore_index=True)
dataframe_to_fasta(df_All_final_ORF, 'Sequence', 'ID',fasta_ORFs_withSP)

#================================ HMMER ANOTATION ===================================================================================================================================================

# We download the pfam database if not present in the folder
if not os.path.exists('./Pfam-A.hmm'):
    print('Download the PFAM database ...')
    urllib.request.urlretrieve("https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz", "Pfam-A.hmm.gz")
    with gzip.open('Pfam-A.hmm.gz', 'rb') as fi:
        with open('Pfam-A.hmm', 'wb') as fo:
            shutil.copyfileobj(fi, fo)

    os.remove("./Pfam-A.hmm.gz")

print('\nRunning HMMER...')

try:
    subprocess.Popen(['hmmsearch','--cut_ga','--cpu','8','--tblout',name+"HMMER_output1_target_table.out","--domtblout",name+"HMMER_output2_domain_table.out","Pfam-A.hmm",name+"_finalORFs.fasta"]).communicate()
except FileNotFoundError:
    print("The hmmsearch program is not installed or was not found by the pipeline")
    sys.exit()


#================================ HMMER FILTER OUTPUT ================================================================================================================================================

def parse_hmmsearch_output(domtblout_filename):
    df_domtblout = pandas.read_csv(domtblout_filename, comment="#", delim_whitespace=True, names=["target name","accession_t","tlen","query name","accession_Q","qlen","E-value","score_1","bias_1","#","of","c-Evalue","i-Evalue","score_2","bias_2","from_1","to_1","from_2","to_2","from_3","to_3","acc","description of target"])
    aggregated_domains = df_domtblout.groupby('target name')['query name'].apply(list).reset_index()
    aggregated_domains['query name'] = aggregated_domains['query name'].apply(lambda x: list(set(x)))
    return aggregated_domains

df_dHMMER = parse_hmmsearch_output((name+'HMMER_output2_domain_table.out'))

#================================ Extract and store information from final output files  ================================================================================================================================================                        

df_FinORFs = fastaToDataframe(name+'_finalORFs.fasta')                      # On recupère les Id:sequence des ORFs ainsi que leur séquence
df_sequenceSignal = fastaToDataframe(name+'_SSeqs.fasta')            # On recupère les id:séquence des signal sequence détecter par signalP
toxDBBlast = pandas.read_csv(name+"_finalORFs.blastsprot.out", sep="\t",names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])                             # Permettra de verifier si un match avec la bdd de toxine a été identifié selon l'ID [0,1,2,3,-2]
finalStat['NB sequence hit with blastp on in-house toxin Database : '] = len(toxDBBlast)                  

#================================ Creation of the final file who merge information about similarity approach and structural approach  ================================================================================================================================================

df_global = pandas.merge(df_FinORFs,df_sequenceSignal,left_on="ID",right_on="ID",how='left')
df_global = pandas.merge(df_global,df_dHMMER[["target name","query name"]],left_on="ID",right_on="target name",how='left')
df_global = pandas.merge(df_global,toxDBBlast[["qseqid","sseqid","pident","length","evalue"]],left_on="ID",right_on="qseqid",how='left')
df_global.rename(columns={'Sequence_x': 'AA_seq',"Sequence_y":"SignalSeq","query name":"HMMERdomain","sseqid":"ID_blast_toxinDB","pident":"PID_blast_toxin","length":"length_blast_toxin","evalue":"Evalue_blast_toxin"}, inplace=True)
df_global = df_global.drop(["target name","qseqid"],axis=1, errors='ignore')
df_global = df_global.drop_duplicates(subset=['ID'], keep='first').reset_index()

df_global['AA_seq'].fillna('')
df_global['SignalSeq'] = df_global['SignalSeq'].fillna('')

def extract_sequence_from_end(long_sequence, short_sequence):
    position = long_sequence.rfind(short_sequence)
    if position == -1:
        return None
    return long_sequence[position + len(short_sequence):]

df_global["MatureSeq"] = df_global.apply(lambda row: extract_sequence_from_end(row['AA_seq'], row['SignalSeq']), axis=1)

df_global.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')

#================================ Repeats Detection   ================================================================================================================================================

def findRepetition(size,seq):
    repetition=[]
    for cdl in range(0,size):
        sub = [seq[i:i+size] for i in range(cdl, len(seq), size)]
        groups = itertools.groupby(sub)
        result = [(label, sum(1 for _ in group)) for label, group in groups]
        for elem, nbRep in result:
            if int(nbRep) >=int(args.RepeatsLength):
               repetition.append((elem,nbRep))
    return repetition

secreted = pandas.read_csv(name+'_secreted_ORFs.alldata',delimiter="\t")

secreted["Repeats1"] = secreted.apply(lambda x: findRepetition(1,x["AA_seq"]),axis=1)
secreted["Repeats2"] = secreted.apply(lambda x: findRepetition(2,x["AA_seq"]),axis=1)
secreted["Repeats3"] = secreted.apply(lambda x: findRepetition(3,x["AA_seq"]),axis=1)
secreted["Repeats"] = secreted["Repeats1"] + secreted["Repeats2"] + secreted["Repeats3"]
secreted['RepeatsTypes'] = secreted['Repeats'].apply(lambda t: [n for (n, _) in t])
secreted['RepeatsLengths'] = secreted['Repeats'].apply(lambda t: [n for (_, n) in t])
secreted['RepeatsLengths'] = [','.join(map(str, l)) for l in secreted['RepeatsLengths']]
secreted['RepeatsTypes'] = [','.join(map(str, l)) for l in secreted['RepeatsTypes']]
secreted = secreted.drop(columns=["Repeats","Repeats1","Repeats2","Repeats3"])
secreted.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')

#================================ Cystein pattern detection ================================================================================================================================================
if args.cysPattern :
    secreted = pandas.read_csv(name+'_secreted_ORFs.alldata',delimiter="\t")

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
        return pattern


    secreted["CysPattern"] = secreted["MatureSeq"].apply(get_cys_pattern)
    secreted.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')

#================================ TPM Quantification Salmon ================================================================================================================================================
if args.Salmon and args.R1ReadsFile :
    secreted = pandas.read_csv(name+'_secreted_ORFs.alldata',delimiter="\t")
    secreted["ID_transcrit"]=secreted.ID.str.split("|",expand=True)[1].str.split(pat=":",n=1,expand=True)[0].str.split(pat="_",n=1,expand=True)[1]
    subprocess.Popen(['salmon','index','-t',transcFilesPath,'-i',('index_'+name)]).communicate()
    subprocess.Popen(['salmon','quant','-i',('index_'+name),'-l','A','-1',args.R1ReadsFile,'-2',args.R2ReadsFile,'-o','output_salmon']).communicate()
    quantiSalmon = pandas.read_csv("./output_salmon/quant.sf",delimiter="\t")
    secreted = pandas.merge(secreted,quantiSalmon[["Name","TPM"]],left_on="ID_transcrit",right_on="Name")
    secreted = secreted.drop(["Name","ID_transcrit_y","ID_transcrit_x"],axis=1, errors='ignore')
    secreted.rename(columns={'tpm': 'TPM_salmon'}, inplace=True)
    secreted.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')

#================================ WoLF PSORT localisation prediction ================================================================================================================================================
if args.wolfPSortPath is not None:
    secreted = pandas.read_csv(name+'_secreted_ORFs.alldata',delimiter="\t")
    awk = "awk '{print $1,$2}'"
    cmd = "{} animal < {}_finalORFs.fasta | {}".format(args.wolfPSortPath,name,awk)

    resultOfWolfPsort = subprocess.Popen(cmd, stdout=subprocess.PIPE,shell=True)
    resultOfWolfPsortToString = StringIO(resultOfWolfPsort.communicate()[0].decode('utf-8'))
    wolfPsortResultInTab = pandas.read_csv(resultOfWolfPsortToString,skiprows=1,delim_whitespace=True,names=["shortID","wolfPSort_localization"])
    wolfPsortResultInTab["ID_transcrit"]=wolfPsortResultInTab.shortID.str.split(pat="|",n=1,expand=True)[1]
    secreted = pandas.merge(secreted,wolfPsortResultInTab[["ID_transcrit","wolfPSort_localization"]],left_on="shortID",right_on="ID_transcrit")
    finalStat["WoLFPSortStat"] = secreted["wolfPSort_localization"].value_counts()
    secreted.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')

#================================ Blast on protein database ================================================================================================================================================
if args.SwissProtDatabase:

    # We download the SwissProt database if not present in the folder
    if not os.path.exists('./uniprot_sprot.fasta'):
        print('Download the SwissProt database ...')
        urllib.request.urlretrieve("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", "uniprot_sprot.fasta.gz")
        with gzip.open('uniprot_sprot.fasta.gz', 'rb') as fi:
            with open('uniprot_sprot.fasta', 'wb') as fo:
                shutil.copyfileobj(fi, fo)

        os.remove("./uniprot_sprot.fasta.gz")


    try:
        CPU = str(mp.cpu_count())
        subprocess.Popen(["makeblastdb","-in","./uniprot_sprot.fasta","-dbtype","prot"]).communicate()
        subprocess.Popen(["blastp","-db","./uniprot_sprot.fasta","-query",name+"_finalORFs.fasta","-out",name+"_blastProtein.blastsprot.out","-outfmt","6","-evalue","1E-5","-max_target_seqs","1","-num_threads",CPU]).communicate()

        secreted = pandas.read_csv(name+'_secreted_ORFs.alldata',delimiter="\t")
        blastResult = pandas.read_csv(name+"_blastProtein.blastsprot.out",delimiter="\t")
        blastResult.columns = ["qseqid","sseqid", "pident", "length","mismatch", "gapopen", "qstart", "qend","sstart", "send", "evalue", "bitscore"]
        #blastResult["qseqid"]=blastResult.qseqid.str.split(pat="|",n=1,expand=True)[1]
        blastResult = blastResult.drop(columns=["length","mismatch","gapopen","qstart","qend","sstart","send","bitscore"])
        secreted = pandas.merge(secreted,blastResult[["qseqid","sseqid","pident","evalue"]],left_on="ID",right_on="qseqid",how="left")
        secreted = secreted.drop_duplicates("ID")
        secreted.rename(columns={'sseqid': 'Protein_blast_ID','pident':'percentage_of_identity_blast_protein','evalue':'Evalue_blast_protein'}, inplace=True)
        secreted = secreted.drop(["qseqid"],axis=1, errors='ignore')
        secreted.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')
        finalStat['NB sequence hit with blastp on protein Database : '] = len(blastResult["qseqid"])
    except FileNotFoundError:
        print("The blast program is not installed or was not found by the pipeline")
        sys.exit()


#================================ Rating proteins ================================================================================================================================================

try:
    secreted = pandas.read_csv(name+'_secreted_ORFs.alldata',delimiter="\t")
    secreted = secreted.assign(Rating='')
    secreted['Rating'] = secreted.apply(lambda row: str(row['Rating'] + 'S') if pandas.notna(row['SignalSeq']) else str(row['Rating'] + '*'), axis=1)
    secreted['Rating'] = secreted.apply(lambda row: str(row['Rating'] + 'B') if (row['ID_blast_toxinDB'] != "nohit" ) else row['Rating'], axis=1)
    if 'CysPattern' in secreted.columns:
        secreted['Rating'] = secreted.apply(lambda row: str(row['Rating'] + 'C') if pandas.notna(row['CysPattern']) else row['Rating'], axis=1)
    if 'TPM' in secreted.columns:
        secreted['Rating'] = secreted.apply(lambda row: str(row['Rating'] + 'T') if (float(row['TPM'])>=float(args.TPMthreshold)) else row['Rating'], axis=1)
    secreted['Rating'] = secreted.apply(lambda row: str(row['Rating'] + 'D') if pandas.notna(row['HMMERdomain']) else row['Rating'], axis=1)
    if 'Protein_blast_ID' in secreted.columns:
        secreted['Rating'] = secreted.apply(lambda row: str(row['Rating'] + '!') if pandas.notna(row['Protein_blast_ID']) and (row['ID_blast_toxinDB'] == "nohit" ) else row['Rating'], axis=1)
    secreted.to_csv(name+'_secreted_ORFs.alldata', index=False,sep='\t')
except:
    print("An error has occurred during sequence rating")
    sys.exit()   



#================================ Backup of pipeline statistics ================================================================================================================================================

def statOfpipeline(outputFilename):
    file = open(outputFilename, 'w')
    for key,value in finalStat.items():
        file.write(key+str(value)+'\n\n')
    file.close()
    
statOfpipeline("stat_Of_pipeline.txt")
