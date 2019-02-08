# ChIP-nexus processing pipeline

## Introduction 

ChIP-nexus pipeline

### Features

* **Flexibility**: Support for `docker`, `singularity` and `Conda`.
* **Portability**: Support for many cloud platforms (Google/DNAnexus) and cluster engines (SLURM/SGE/PBS).
* **Resumability**: [Resume](utils/qc_jsons_to_tsv/README.md) a failed workflow from where it left off.
* **User-friendly HTML report**: tabulated quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures).
* **Genomes**: Pre-built database for GRCh38, hg19, mm10, mm9 and additional support for custom genomes.

## Installation and tutorial

This pipeline supports many cloud platforms and cluster engines. It also supports `docker`, `singularity` and `Conda` to resolve complicated software dependencies for the pipeline. A tutorial-based instruction for each platform will be helpful to understand how to run pipelines. There are special instructions for two major Stanford HPC servers (SCG4 and Sherlock).

* Cloud platforms
  * Web interface
    * [DNAnexus Platform](docs/tutorial_dx_web.md)
  * CLI (command line interface)
    * [Google Cloud Platform](docs/tutorial_google.md)
    * [DNAnexus Platform](docs/tutorial_dx_cli.md)
* Stanford HPC servers (CLI)
  * [Stanford SCG4](docs/tutorial_scg.md)
  * [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)
* Cluster engines (CLI)
  * [SLURM](docs/tutorial_slurm.md)
  * [Sun GridEngine (SGE/PBS)](docs/tutorial_sge.md)
* Local computers (CLI)
  * [Local system with `singularity`](docs/tutorial_local_singularity.md)
  * [Local system with `docker`](docs/tutorial_local_docker.md)
  * [Local system with `Conda`](docs/tutorial_local_conda.md)

## Input JSON file

[Input JSON file specification](docs/input.md)

## Output directories

[Output directory specification](docs/output.md)

## Useful tools

There are some useful tools to post-process outputs of the pipeline.

### qc_jsons_to_tsv

[This tool](utils/qc_jsons_to_tsv/README.md) recursively finds and parses all `qc.json` (pipeline's [final output](docs/example_output/v1.1.5/qc.json)) found from a specified root directory. It generates a TSV file that has all quality metrics tabulated in rows for each experiment and replicate. This tool also estimates overall quality of a sample by [a criteria definition JSON file](utils/qc_jsons_to_tsv/criteria.default.json) which can be a good guideline for QC'ing experiments.

### resumer

[This tool](utils/resumer/README.md) parses a metadata JSON file from a previous failed workflow and generates a new input JSON file to start a pipeline from where it left off.

### ENCODE downloader

[This tool](https://github.com/kundajelab/ENCODE_downloader) downloads any type (FASTQ, BAM, PEAK, ...) of data from the ENCODE portal. It also generates a metadata JSON file per experiment which will be very useful to make an input JSON file for the pipeline.