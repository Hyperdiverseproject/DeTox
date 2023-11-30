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
