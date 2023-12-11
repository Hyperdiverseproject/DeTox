# DeTox pipeline

## Installation


### 1. Clone Repository
- Clone the repository to your local machine.
  ```bash
  git clone <repository_url>
  ```

### 2. Create and Activate Conda Environment
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
  mamba activate DeToX && phobius-register phobius101_linux.tar.gz
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
```bash
snakemake --snakefile /path/to/snakefile -r all -j 8 --configfile config.yaml
```

### Configuration of settings and options
To use this project, you need to fill the `config.yaml` file with the following parameters:

- `transcriptome`: the path to the transcriptome file in FASTA format. If this parameter is not provided, the assembly will be performed using the `R1` and `R2` parameters.
- `basename`: the basename for the output files.
- `memory`: the memory in GB to use for the assembly. This parameter might change later.
- `threads`: the number of logic threads, not physical cores, to use for the assembly. For clarity, this parameter should be equal to or less than the number of cores available on your machine.
- `R1`: the path to the forward reads file in FASTQ format. This parameter is required if `transcriptome` is not provided.
- `R2`: the path to the reverse reads file in FASTQ format. This parameter is required if `transcriptome` is not provided.
- `adapters`: the path to the adapters file in FASTA format. This parameter is required if `transcriptome` is not provided.
- `contaminants`: the path to the contaminants file in FASTA format. This parameter is mandatory.
- `toxin_db`: the path to the toxin database file in FASTA format. This parameter is mandatory.
- `toxins_evalue`: the e-value threshold for toxin annotation using BLAST.
- `contamination_evalue`: the e-value threshold for contamination removal using BLAST.
- `clustering_threshold`: the clustering threshold for CD-HIT-EST.
- `maxlen`: the maximum length of transcripts to keep.
- `signalp_dvalue`: the D-value threshold for SignalP.
- `signalp_path`: the path to the SignalP executable.
- `wolfPsort_path`: the path to the WoLF PSORT executable (WoLFPSort/bin/runWolfPsortSummary). This parameter is required if `wolfpsort` is set to `True`.
- `quant`: a boolean option to perform transcript quantification using Salmon. If set to `True`, the `R1` and `R2` parameters are required.
- `TPMthreshold`: The TPM threshold for a sequence to be flagged "T".
- `wolfpsort`: a boolean option to perform subcellular localization prediction using WoLF PSORT. If set to `True`, the `wolfPsort_path` parameter is required.
- `swissprot`: a boolean option to perform functional annotation using SwissProt. If set to `True`, the `swissprot_evalue` parameter is required.
- `swissprot_evalue`: the e-value threshold for SwissProt annotation using BLAST.
- `cys_pattern`: a boolean option to perform cysteine pattern analysis.

