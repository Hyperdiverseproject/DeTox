# DeTox pipeline

## install
+ clone this repository
+ switch to snakemake_port branch
  + `git switch snakemake_port` or `git checkout snakemake_port`
+ create the conda environment (here i use mamba for faster setup)
  + `mamba env create -f toxo_env.yml && mamba activate DeToX`
+ install licensed software
  + **signalp**
    + obtain signalP5 from [the official website](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=5.0&packageversion=5.0b&platform=Linux)
    + follow the provided readme file for system-wide installation
  + **phobius**
    + obtain a copy of phobius from [here](https://phobius.sbc.su.se/data.html)
    + the downloaded file should be named `phobius101_linux.tar.gz`, move to its directory within a terminal
    + the environment includes a practical script for installing phobius: 
      + `mamba activate DeToX && phobius-register phobius101_linux.tar.gz`
  + **wolfPSoRT** (TODO: automatic install this)
    +  clone [this repo](https://github.com/fmaguire/WoLFPSort)
    +  copy the path to `runWolfPsortSummary` (e.g., `~/WoLFPSort/bin/runWolfPsortSummary`) in the proper section of `config.yaml`

## usage
+ activate the environment
  + `mamba activate DeToX`
  + copy `config.yaml` to your working directory
  + populate it 
  + run the pipeline
    + `snakemake --snakefile /path/to/snakefile -r all -j 8 --configfile config.yaml` 

## snakemake port

## stuff to be done:

- [x] port C pattern recognition
  - [x] split mature peptides
- [x] port repetitive patterns flagging
- [x] port merging of all results in a single table
- [ ] port scoring

- [ ] automatic install of wolfpsort
- [ ] find a wolfpsort substitute
- [?] create a wrapper for the script

## current status of the port graph:
![DAG](dag.svg)
