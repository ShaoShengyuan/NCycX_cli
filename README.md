# NCycX: Nitrogen cycling genes expanded database and high-resolution explorer for microbiomes

[![Web Server](https://img.shields.io/badge/Web_Server-ncycx.com-blue?style=flat&logo=google-chrome)](https://ncycx.com)

Welcome to the **NCycX** Command-Line Interface (CLI) repository! 

Microbial nitrogen cycling drives global biogeochemical processes and plays a crucial role in maintaining ecosystem functioning and agricultural sustainability. However, comprehensive reference resources for microbial nitrogen-cycling genes remain limited. Therefore, we established **NCycX** (*Nitrogen Cycling genes eXpanded database and high-resolution eXplorer for microbiomes*), a comprehensive and advanced resource for the systematic profiling and exploration of nitrogen-cycling and transport genes.

The current release includes a large-scale reference dataset, comprising 153 core functional genes (encompassing inorganic/organic nitrogen metabolism and transporters) derived from more than 110,000 microbial genomes, complemented by protein domain profiles and partial 3D structural data. 

**Web Server**: For user-friendly web interfaces, sequence annotation, and homology-based searches, please visit our official web portal at **[https://ncycx.com](https://ncycx.com)**. 

**CLI Tool**: This repository provides the downloadable datasets and local command-line tools for offline, large-scale, and high-throughput data processing.

---

## Installation & Dependencies

We recommend using `conda` to manage the environment and install the required dependencies.

```bash
# 1. Create and activate a new conda environment
conda create -n ncycx -y
conda activate ncycx

# 2. Install alignment and sequence analysis tools
conda install -c conda-forge -c bioconda diamond blast -y
conda install -c conda-forge -c bioconda hmmer -y
```

> **Important Note on Database Files:**
> Due to file size limits on GitHub, the core database files (e.g., `.faa`, `.dmnd`) in the `data/` directory may be compressed (e.g., `.gz` format) or hosted externally (like Zenodo). 
> **Before running the pipeline, please ensure you have decompressed the database files:**
> ```bash
> cd data/
> gunzip *.gz
> cd ..
> ```

---

## Quick Start (Examples)

### 1. Configuration
Before running the pipeline, you need to set up your configuration file (`.ctl`). Open the `.ctl` file (e.g., `ncycx_core.ctl`) and modify the following key sections to match your local environment:
* `# Main Workflow Parameters`: Set sequence format, input directory, etc.
* `# Database and Mapping File Paths`: Ensure these point to the correct files in your `data/` directory.
* `# Executable Paths`: Set the absolute paths for external software (Diamond, BLAST, HMMER) if they are not already in your system's `$PATH`.

### 2. Generate Sample Information File (SI.txt)
If you are analyzing multiple `.faa` files, you can quickly generate a sample information file using `seqkit`:

```bash
cd examples/

# Install seqkit if you haven't already
conda install -c conda-forge -c bioconda seqkit -y

# Generate SI.txt based on .faa files in the current directory
seqkit stat -T *.faa | tail -n +2 | sed 's/\.faa//g' | cut -f 1,4 | xargs -I{} basename {} > SI.txt
```

### 3. Run the CLI Pipeline
Execute the Python wrapper script, specifying your configuration file and desired output path:

```bash
# Ensure you specify the correct paths for your input directory and config file
python NCycX_cli.py -c ncycx_core.ctl -o prot_diamond.out
```

---

## Command-Line Arguments (Help)

The `NCycX_cli.py` script supports an initial screening + HMM fine-filtering pipeline. You can use command-line arguments to easily override the parameters defined in your `.ctl` configuration file.

```text
usage: python NCycX_cli.py [-h] -c CONFIG [-d DIR] [-m METHOD] [-f FILETYPE] [-s SEQTYPE]
                           [-si SAMPLEINFO] [-rs RANDOMSAMPLING] [-o OUTFILE] [--hmm-filter]

NCycX DB wrapper: Supports initial screening + HMM fine-filtering pipeline

optional arguments:
  -h, --help            Show this help message and exit
  -c CONFIG, --config CONFIG
                        Path to the configuration file [Required] (e.g., config.ctl)
  -d DIR, --dir DIR     Override: Directory containing sequence files
  -m METHOD, --method METHOD
                        Override: Alignment software (diamond/usearch/blast)
  -f FILETYPE, --filetype FILETYPE
                        Override: Sequence file extension (e.g., fasta, faa)
  -s SEQTYPE, --seqtype SEQTYPE
                        Override: Sequence type (prot or nucl)
  -si SAMPLEINFO, --sampleinfo SAMPLEINFO
                        Override: Tab-separated file containing sample info
  -rs RANDOMSAMPLING, --randomsampling RANDOMSAMPLING
                        Override: Rarefaction size (set to 0 to skip)
  -o OUTFILE, --outfile OUTFILE
                        Override: Output path for the final abundance table
  --hmm-filter          Whether to perform HMM domain filtering after alignment
```