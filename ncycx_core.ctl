# ==========================================
# Main Workflow Parameters (核心工作流参数)
# ==========================================

# Directory containing your sequence files (e.g., ./data)
# 指定存放输入序列文件的文件夹路径
WORKDIR=/home/shaosy/264NDB/NCycX/examples/faa

# Path to the tab-delimited file containing sample names and sequence numbers
# 包含样本名和序列数量的制表符分隔文件
SAMPLEINFO=/home/shaosy/264NDB/NCycX/examples/SI_prot.txt

# Path to the output functional profile result file
# 最终生成的丰度结果表路径
OUTFILE=./abundance.out

# Alignment method to use. Options: diamond, usearch, blast. (Default: diamond)
# 使用的比对软件，默认 diamond
METHOD=diamond

# File extension of your sequence files. 
# Options: fastq, fastq.gz, fasta, fasta.gz, fq, fq.gz, fa, fa.gz, faa, faa.gz
# 输入文件的后缀类型
FILETYPE=faa

# Sequence type of your input files. Options: prot (protein), nucl (nucleotide)
# 输入序列的类型：prot(蛋白质) 或 nucl(核酸)
SEQTYPE=prot

# Random sampling size (integer). 
# 设为 0，表示跳过抽平，输出最真实的绝对序列计数 (Absolute Counts)
# 如果留空，则依旧会自动寻找最小样本数进行抽平
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
