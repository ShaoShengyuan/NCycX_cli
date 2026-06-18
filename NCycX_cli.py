#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import argparse
import pandas as pd
from collections import defaultdict

def parse_config(ctl_path):
    """解析 .ctl 配置文件"""
    config = {}
    with open(ctl_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, val = line.split('=', 1)
                key = key.strip()
                val = val.strip().replace('~', os.path.expanduser('~'))
                if val:
                    config[key] = val
    return config

def run_perl_script(args):
    """调用底层的 Perl 脚本执行初步比对和定量"""
    print("="*50)
    print("🚀 [Step 1] Running Perl script for initial alignment...")
    print("="*50)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    perl_script_path = os.path.join(script_dir, "Ngenes_2606.PL")
    
    if not os.path.exists(perl_script_path):
        sys.exit(f"\n[Error] Underlying Perl script not found! Please ensure {perl_script_path} exists.")

    cmd = ["perl", perl_script_path, "-c", args.config]
    
    if args.dir:             cmd.extend(["-d", args.dir])
    if args.method:          cmd.extend(["-m", args.method])
    if args.filetype:        cmd.extend(["-f", args.filetype])
    if args.seqtype:         cmd.extend(["-s", args.seqtype])
    if args.sampleinfo:      cmd.extend(["-si", args.sampleinfo])
    if args.randomsampling is not None: cmd.extend(["-rs", args.randomsampling])
    if args.outfile:         cmd.extend(["-o", args.outfile])
        
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit("\n[Error] Perl script execution failed. Pipeline terminated.")
    print("\n✅ Initial alignment by Perl script completed!")

def extract_fasta_for_hmm(config):
    """提取 FASTA 序列，并记录所有初步被注释的序列映射"""
    workdir = config.get('WORKDIR')
    filetype = config.get('FILETYPE')
    method = config.get('METHOD', 'diamond')
    map_file = config.get('MAP_FILE')
    
    print("\n" + "="*50)
    print("🚀 [Step 2] Extracting FASTA sequences based on initial alignment results...")
    print("="*50)

    id2gene = {}
    with open(map_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                id2gene[parts[0]] = parts[1]

    hmm_fasta_dir = os.path.join(workdir, "hmm_temp_fasta")
    os.makedirs(hmm_fasta_dir, exist_ok=True)
    
    sample_files = glob.glob(os.path.join(workdir, f"*{filetype}"))
    if not sample_files:
        sys.exit(f"[Error] No sequencing files with format {filetype} found in {workdir}.")
    if filetype.endswith('.gz') or filetype.startswith('fastq') or filetype.startswith('fq'):
        sys.exit("\n[Error] The --hmm-filter pipeline currently ONLY supports uncompressed FASTA files (.fasta, .faa, .fa). Please decompress or convert your input files before using the HMM filter.")

    gene_file_handles = {}
    extracted_count = 0
    
    # 【新增】：用于追踪所有被初步注释的序列
    # 格式： seq_to_gene["sample_name::query_id"] = "gene_type"
    seq_to_gene = {} 
    
    for faa_path in sample_files:
        filename = os.path.basename(faa_path)
        sample_name = filename.replace(f".{filetype}", "")
        aln_path = os.path.join(workdir, f"{sample_name}.{method}")
        
        if not os.path.exists(aln_path):
            continue
            
        query2gene = {}
        with open(aln_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                query_id, db_id = parts[0], parts[1]
                gene_type = id2gene.get(db_id)
                if gene_type and query_id not in query2gene:
                    query2gene[query_id] = gene_type

        with open(faa_path, 'r') as f:
            write_flag = False
            current_gene = None
            for line in f:
                if line.startswith('>'):
                    query_id = line.strip()[1:].split()[0]
                    if query_id in query2gene:
                        write_flag = True
                        current_gene = query2gene[query_id]
                        full_seq_name = f"{sample_name}::{query_id}"
                        
                        # 记录到全局字典中，方便最后对比
                        seq_to_gene[full_seq_name] = current_gene
                        
                        if current_gene not in gene_file_handles:
                            out_path = os.path.join(hmm_fasta_dir, f"{current_gene}.fasta")
                            gene_file_handles[current_gene] = open(out_path, 'w')
                        
                        gene_file_handles[current_gene].write(f">{full_seq_name}\n")
                        extracted_count += 1
                    else:
                        write_flag = False
                elif write_flag:
                    gene_file_handles[current_gene].write(line)

    for fh in gene_file_handles.values():
        fh.close()
        
    print(f"✅ Successfully extracted {extracted_count} sequences and sorted into: {hmm_fasta_dir}")
    return hmm_fasta_dir, seq_to_gene

def run_hmmscan(config, hmm_fasta_dir):
    """对提取的 FASTA 运行 hmmscan"""
    workdir = config.get('WORKDIR')
    hmmscan = config.get('HMMSCAN_PATH')
    hmm_db = config.get('HMM_DB')
    threads = config.get('HMM_THREADS', '5')
    
    hmm_result_dir = os.path.join(workdir, "hmmscan_result")
    os.makedirs(hmm_result_dir, exist_ok=True)
    
    print("\n" + "="*50)
    print("🚀 [Step 3] Executing hmmscan...")
    print("="*50)
    
    fasta_files = glob.glob(os.path.join(hmm_fasta_dir, "*.fasta"))
    total = len(fasta_files)
    
    for i, fasta_path in enumerate(fasta_files, 1):
        gene_type = os.path.basename(fasta_path).replace(".fasta", "")
        print(f"[{i}/{total}] Scanning {gene_type}...")
        
        tbl_out = os.path.join(hmm_result_dir, f"{gene_type}.tbl")
        dom_out = os.path.join(hmm_result_dir, f"{gene_type}.dom")
        
        cmd = [
            hmmscan, "--cpu", str(threads),
            "--tblout", tbl_out,
            "--domtblout", dom_out,
            "-o", os.path.join(hmm_result_dir, f"{gene_type}.log"),
            hmm_db, fasta_path
        ]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f"✅ hmmscan execution completed! Results saved to: {hmm_result_dir}")
    return hmm_result_dir

def hmm_filter_and_recalc(config, hmm_result_dir, seq_to_gene):
    """核心：使用你提供的可靠 Pandas 逻辑过滤，并输出丰度表与追踪表"""
    tsv_path = config.get('HIGH_FREQ_TSV')
    original_outfile = config.get('OUTFILE', './abundance.txt')
    
    if original_outfile.endswith('.txt'):
        final_outfile = original_outfile.replace('.txt', '_hmm_filtered.txt')
        validation_file = original_outfile.replace('.txt', '_hmm_validation.tsv')
    else:
        final_outfile = original_outfile + '_hmm_filtered.txt'
        validation_file = original_outfile + '_hmm_validation.tsv'
    
    print("\n" + "="*50)
    print("🚀 [Step 4] Fine-filtering and recalculating abundance based on high_freq_domain.tsv...")
    print("="*50)
    
    # 1. 采用你的逻辑加载字典
    print("Loading reference Domain dictionary...")
    all_high_freq = pd.read_csv(tsv_path, sep='\t')
    gene_to_accessions = dict()
    for gene_type, group in all_high_freq.groupby('gene_type'):
        gene_to_accessions[gene_type] = group['accession'].tolist()
        
    filtered_counts = defaultdict(lambda: defaultdict(int))
    all_samples = set()
    
    # 用于记录哪些序列成功活过了 HMM 筛选
    passed_sequences = set()
    
    tbl_files = glob.glob(os.path.join(hmm_result_dir, "*.tbl"))
    
    for tbl_path in tbl_files:
        tbl_file = os.path.basename(tbl_path)
        gene_type = os.path.splitext(tbl_file)[0]
        
        if gene_type not in gene_to_accessions:
            print(f"❌ {gene_type} not found in all_high_freq.tsv, skipping...")
            continue
            
        reference_accessions = gene_to_accessions[gene_type]
        
        # 3. 完全复用你的 .tbl 读取逻辑
        try:
            domain_tbl = pd.read_csv(
                tbl_path, header=None, skiprows=3, skipfooter=10,
                delimiter=r'\s+', engine='python', usecols=range(17)
            )
        except Exception as e:
            # 如果文件为空或者读取失败，安静跳过
            continue
        # 脱除版本号 (例如 PF00001.21 变成 PF00001)，增强匹配鲁棒性
        col1_clean = domain_tbl[1].astype(str).str.split('.').str[0]
            
        # 4. 核心过滤：基于清理后的 Accession 列进行精准匹配 
        filtered_rows = domain_tbl[col1_clean.isin(reference_accessions)]
        unique_sequences = set(filtered_rows[2].tolist())
        
        # 将通过筛选的序列存入总集合
        passed_sequences.update(unique_sequences)
        
        for seq_name in unique_sequences:
            if "::" in seq_name:
                sample_name, _ = seq_name.split("::", 1)
                filtered_counts[sample_name][gene_type] += 1
                all_samples.add(sample_name)
    
    print(f"\n📊 Filtering statistics: Initial screen {len(seq_to_gene)} seqs -> HMM retained {len(passed_sequences)} seqs.")
    
    # ================= 重新生成丰度矩阵 =================
    all_genes = set()
    for sample in filtered_counts:
        all_genes.update(filtered_counts[sample].keys())
        
    sorted_samples = sorted(list(all_samples))
    sorted_genes = sorted(list(all_genes))
    
    with open(final_outfile, 'w') as out:
        out.write("# HMM Filtered Profile\n")
        out.write("Gene\t" + "\t".join(sorted_samples) + "\n")
        
        for gene in sorted_genes:
            row = [gene]
            for sample in sorted_samples:
                count = filtered_counts.get(sample, {}).get(gene, 0)
                row.append(str(count))
            out.write("\t".join(row) + "\n")

    # ================= 生成验证/追踪表格 =================
    print(f"📝 Generating sequence tracking validation table: {validation_file}...")
    with open(validation_file, 'w') as vf:
        vf.write("Sequence_ID\tAnnotated_Gene\tHMM_Passed\n")
        for full_seq_name, gene in seq_to_gene.items():
            status = "Yes" if full_seq_name in passed_sequences else "No"
            vf.write(f"{full_seq_name}\t{gene}\t{status}\n")

    print(f"✅ HMM fine-filtering completed successfully!")
    print(f"   - Final abundance table: {final_outfile}")
    print(f"   - Sequence validation table: {validation_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NCycX DB wrapper: Supports initial screening + HMM fine-filtering pipeline")
    
    parser.add_argument("-c", "--config", required=True, help="Path to the configuration file [Required] (e.g., config.ctl)")
    parser.add_argument("-d", "--dir", help="Override: Directory containing sequence files")
    parser.add_argument("-m", "--method", help="Override: Alignment software (diamond/usearch/blast)")
    parser.add_argument("-f", "--filetype", help="Override: Sequence file extension (e.g., fasta, faa)")
    parser.add_argument("-s", "--seqtype", help="Override: Sequence type (prot or nucl)")
    parser.add_argument("-si", "--sampleinfo", help="Override: Tab-separated file containing sample info")
    parser.add_argument("-rs", "--randomsampling", help="Override: Rarefaction size (set to 0 to skip)")
    parser.add_argument("-o", "--outfile", help="Override: Output path for the final abundance table")
    parser.add_argument("--hmm-filter", action="store_true", help="Whether to perform HMM domain filtering after alignment")
    
    args = parser.parse_args()
    
    config = parse_config(args.config)
    
    if args.dir: config['WORKDIR'] = args.dir
    if args.method: config['METHOD'] = args.method
    if args.filetype: config['FILETYPE'] = args.filetype
    if args.seqtype: config['SEQTYPE'] = args.seqtype
    if args.sampleinfo: config['SAMPLEINFO'] = args.sampleinfo
    if args.randomsampling is not None: config['RANDOMSAMPLING'] = args.randomsampling
    if args.outfile: config['OUTFILE'] = args.outfile
    
    run_perl_script(args)
    
    if args.hmm_filter:
        hmm_fasta_dir, seq_to_gene = extract_fasta_for_hmm(config)
        hmm_result_dir = run_hmmscan(config, hmm_fasta_dir)
        hmm_filter_and_recalc(config, hmm_result_dir, seq_to_gene)
    else:
        print(f"\n🎉 Analysis completed! Initial results saved to: {config.get('OUTFILE', './abundance.txt')}")
