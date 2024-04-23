# DeTox pipeline

## Installation


### 1. Clone Repository
- Clone the repository to your local machine.
  ```bash
  git clone <repository_url>
  ```

### 2. Create and Activate Conda Environment

#### 2.1 Installation of Conda

- Before you begin, make sure you have a recent version of Python installed on your system. Conda is typically distributed with the Anaconda or Miniconda distribution.

1. **Download the Anaconda Distribution:**
   - Visit the official Anaconda website at [https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
   - Select the appropriate version for your operating system (Windows, macOS, Linux) and download the installation file.

2. **Install Anaconda or Miniconda:**
   - Follow the installation instructions provided on the official Conda page.
   - After installation, you should have Conda installed on your system.

#### 2.2 Installation of Mamba

Mamba is an alternative to Conda that offers improved performance for environment and package management.

1. **Install Mamba:**
   - After installing Conda, open your terminal or command prompt.
   - Type the following command to install Mamba:
     ```
     conda install -c conda-forge mamba
     ```

2. **Verify the Installation:**
   - To ensure that Mamba has been successfully installed, type the following command:
     ```
     mamba --version
     ```
   - You should see the Mamba version displayed, confirming a successful installation.

#### 2.3 Create and Activate Conda Environment

- Utilize Mamba for efficient environment setup and activation. 
  ```bash
  mamba env create -f toxo_env.yml && mamba activate DeToX
  ```

### 3. Install Licensed Software

#### SignalP
- Download SignalP 5.0 from the [official website](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=5.0&packageversion=5.0b&platform=Linux).
- Follow the instructions in the provided README for system-wide installation.

#### Phobius
- Download Phobius from [this link](https://phobius.sbc.su.se/data.html).
- Navigate to the download directory and register Phobius using the provided script:
  ```bash
  mamba activate DeToX && phobius-register phobius101_linux.tgz
  ```

#### WoLF PSORT (Pending Automatic Installation)
- Clone the [WoLF PSORT repository](https://github.com/fmaguire/WoLFPSort).
- Specify the path to `runWolfPsortSummary` in the `config.yaml` file.

---

## Usage

### Activate Environment and Prepare Configuration
1. Activate the `DeToX` environment.
   ```bash
   mamba activate DeToX
   ```
2. Copy `config.yaml` to your working directory and populate it with required information.

### Run Pipeline
Execute the pipeline using Snakemake.

- This Snakemake command utilizes the -j option to specify the maximum number of threads (parallel jobs) for rule execution. The specific number, in this case, is set to 8 (-j 8). Additionally, the number of threads can be influenced by the "threads" parameter in the config.yaml file. Ensure the desired thread count is configured in the "threads" parameter of the config.yaml file before executing the command.
```bash
snakemake --snakefile /path/to/snakefile -r all -j 8 --configfile config.yaml
```

### Configuration of settings and options
To use this project, you need to fill the `config.yaml` file with the following parameters:

- `transcriptome`: the path to the transcriptome file in FASTA format. If this parameter is not provided, the assembly will be performed using the `R1` and `R2` parameters.
- `output_dir`: absolute path of the directory used for output files. If no path is provided the files will be created in the current directory.
- `basename`: the basename for the output files.
- `memory`: the memory in GB to use for the assembly.
- `threads`: the number of logic threads, not physical cores, to use for the assembly. For clarity, this parameter should be equal to or less than the number of cores available on your machine.
- `R1`: the path to the forward paired-end or single-end reads file in FASTQ format. This parameter is required if `transcriptome` is not provided.
- `R2`: the path to the reverse paired-end reads file in FASTQ format. This parameter is required if `transcriptome` is not provided and data are paired-end.
- `adapters`: the path to the adapters file in FASTA format. This parameter is required if `transcriptome` is not provided are paired end.
- `contaminants`: the path to the contaminants file in FASTA format. This parameter is mandatory.
- `toxin_db`: the path to the toxin database file in FASTA format. This parameter is mandatory.
- `pfam_db_path`: the path to the Pfam database file in .hmm format. If no path is provided, the database is downloaded automatically.
- `swissprot_db_path`: the path to the SwissProt database file in fasta.gz format. If no path is provided, the database is downloaded automatically.
- `toxins_evalue`: the e-value threshold for toxin annotation using BLAST.
- `contamination_evalue`: the e-value threshold for contamination removal using BLAST.
- `clustering_threshold`: the clustering threshold for CD-HIT.
- `maxlen`: the maximum length of transcripts to keep.
- `signalp_dvalue`: the D-value threshold for SignalP.
- `wolfPsort_path`: the path to the WoLF PSORT executable (WoLFPSort/bin/runWolfPsortSummary). This parameter is required if `wolfpsort` is set to `True`.
- `quant`: a boolean option to perform transcript quantification using Salmon. If set to `True`, the `R1` and `R2` parameters are required.
- `TPMthreshold`: The TPM threshold for a sequence to be flagged "T".
- `wolfpsort`: a boolean option to perform subcellular localization prediction using WoLF PSORT. If set to `True`, the `wolfPsort_path` parameter is required.
- `swissprot`: a boolean option to perform functional annotation using SwissProt. If set to `True`, the `swissprot_evalue` parameter is required.
- `swissprot_evalue`: the e-value threshold for SwissProt annotation using BLAST.
- `cys_pattern`: a boolean option to perform cysteine pattern analysis.

## Citation

If you use DeTox in your research, please cite our article:

Allan Ringeval, Sarah Farhat, Alexander Fedosov, Marco Gerdol, Samuele Greco, Lou Mary, Maria Vittoria Modica, Nicolas Puillandre, DeTox: a pipeline for the detection of toxins in venomous organisms, Briefings in Bioinformatics, Volume 25, Issue 2, March 2024, bbae094, https://doi.org/10.1093/bib/bbae094
