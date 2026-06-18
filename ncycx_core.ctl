# ==========================================
# Main Workflow Parameters (核心工作流参数)
# ==========================================

# Directory containing your sequence files (e.g., ./data)
WORKDIR=/home/shaosy/264NDB/NCycX/examples/faa

# Path to the tab-delimited file containing sample names and sequence numbers
SAMPLEINFO=/home/shaosy/264NDB/NCycX/examples/SI_prot.txt

# Path to the output functional profile result file
OUTFILE=./abundance.out

# Alignment method to use. Options: diamond, usearch, blast. (Default: diamond)
METHOD=diamond

# File extension of your sequence files. 
# Options: fastq, fastq.gz, fasta, fasta.gz, fq, fq.gz, fa, fa.gz, faa, faa.gz
FILETYPE=faa

# Sequence type of your input files. Options: prot (protein), nucl (nucleotide)
SEQTYPE=prot

# Random sampling size (integer). 
# Set to 0 to skip smoothing and output the most accurate absolute sequence counts (Absolute Counts)
# If left blank, the system will still automatically determine the minimum sample size for smoothing 
RANDOMSAMPLING=0

# ==========================================
# Database and Mapping File Paths (数据库文件路径)
# ==========================================

# Path to the amino acid database file (.faa) (core)
DB_FAA=/home/shaosy/264NDB/NCycX/data/ncycx_2604_core.faa

# Prefix for the Diamond database (Diamond automatically adds .dmnd during makedb)
DB_PREFIX=/home/shaosy/264NDB/NCycX/data/ncycx_2604_core

# Path to the gene mapping file (.map)
MAP_FILE=/home/shaosy/264NDB/NCycX/data/ncycx_2604_core.map

# ==========================================
# Executable Paths (软件执行路径)
# ==========================================
# Use absolute paths or ~/ for the home directory

DIAMOND_PATH=~/miniconda3/envs/ngenesdb/bin/diamond
USEARCH_PATH=/home/shaosy/biosoft/usearch/usearch11.0.667_i86linux64
MAKEBLASTDB_PATH=~/miniconda3/envs/ngenesdb/bin/makeblastdb
BLASTP_PATH=~/miniconda3/envs/ngenesdb/bin/blastp
BLASTX_PATH=~/miniconda3/envs/ngenesdb/bin/blastx

# ==========================================
# Tool-Specific Thresholds & Parameters (比对阈值参数)
# ==========================================

# Diamond settings
DIAMOND_EVALUE=0.0001
DIAMOND_K=1
DIAMOND_THREADS=5

# Usearch settings
USEARCH_ID=0.3

# BLAST+ settings
BLAST_EVALUE=1e-4
BLAST_MAX_TARGET_SEQS=1

# ==========================================
# HMM Filter Settings (HMM 精筛参数)
# ==========================================
# hmmscan executable path
HMMSCAN_PATH=~/miniconda3/envs/ngenesdb/bin/hmmscan

# Reference HMM database path (must be hmmpressed)
HMM_DB=/home/shaosy/236CuMMOs/hmmscan/Pfam-A.hmm

# High frequency domain reference TSV path
HIGH_FREQ_TSV=/home/shaosy/264NDB/NCycX/data/reference_domain.tsv

# Number of CPU threads for hmmscan
HMM_THREADS=10
