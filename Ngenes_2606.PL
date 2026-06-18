#!/usr/bin/perl
use strict;
use List::Util qw(shuffle);
use Getopt::Long;

my ( $workdir, $method, $outfile, $seqtype, $filetype, $sampleinfo, $randomsampling, $config_file );

# 命令行参数现在作为可选的“覆盖(Override)”项
GetOptions(
  "c=s"  => \$config_file,    ## configuration file (Required)
  "d=s"  => \$workdir,        ## override: sequence file directory
  "m=s"  => \$method,         ## override: method (diamond/usearch/blast)
  "f=s"  => \$filetype,       ## override: file type
  "s=s"  => \$seqtype,        ## override: prot/nucl
  "si=s" => \$sampleinfo,     ## override: sample info file
  "rs=s" => \$randomsampling, ## override: random sampling size
  "o=s"  => \$outfile         ## override: output file
);

# 必须提供配置文件
if (!defined $config_file || !-e $config_file) {
  &PrintHelp();
  die "\n[Error] Configuration file is required. Please use: perl Ngenes_2604_full.PL -c config.ctl\n";
}

## 1. 解析配置文件 (Parse configuration file)
my %config;
open(CFG, "<", $config_file) or die "Cannot open configuration file: $config_file\n";
while (<CFG>) {
  chomp;
  s/^\s+//; s/\s+$//; # 去除首尾空格 (Remove leading/trailing spaces)
  next if /^#/ || /^$/; # 跳过注释和空行 (Skip comments and empty lines)
  if (/^([^=]+)=(.*)$/) {
    my ($key, $value) = ($1, $2);
    $key =~ s/^\s+|\s+$//g;
    $value =~ s/^\s+|\s+$//g;
    $value =~ s/^~/$ENV{"HOME"}/; # 展开波浪号为家目录 (Expand tilde to home dir)
    $config{$key} = $value if $value ne ""; # 只有当值不为空时才存入哈希 (Store only if not empty)
  }
}
close CFG;

## 2. 合并参数 (Merge CLI arguments with Config parameters)
# 如果命令行提供了参数，则优先使用命令行；否则使用配置文件中的值
$workdir        = $workdir        || $config{WORKDIR};
$method         = $method         || $config{METHOD} || "diamond"; 
$filetype       = $filetype       || $config{FILETYPE};
$seqtype        = $seqtype        || $config{SEQTYPE};
$sampleinfo     = $sampleinfo     || $config{SAMPLEINFO};
$randomsampling = $randomsampling || $config{RANDOMSAMPLING};
$outfile        = $outfile        || $config{OUTFILE};

## 3. 参数基础校验 (Basic parameter validation)
if ( !defined $workdir
  || !defined $filetype
  || !defined $outfile
  || !defined $seqtype
  || !defined $sampleinfo
  || $method !~ /^diamond|usearch|blast$/
  || $filetype !~ /^fastq|fastq.gz|fasta|fasta.gz|fq|fq.gz|faa|faa.gz|fa|fa.gz|fna|fna\.gz|ffn|ffn\.gz$/ )
{
  &PrintHelp();
  die "\n[Error] Missing or incorrect core parameters. Please check your config file.\n";
}

$workdir =~ s/\/$//; # 去除路径末尾的斜杠 (Remove trailing slash)

## 4. 开始比对分析 (Start alignment analysis)
my (@files);
if ( $method eq "diamond" ) {
  my $diamond = $config{DIAMOND_PATH};
  die "Please specify the location of diamond in your config!" if !-e $diamond;
  @files = glob("$workdir/*$filetype");
  
  my $db_prefix  = $config{DB_PREFIX};
  my $db_faa     = $config{DB_FAA};
  my $diamond_db = "$diamond makedb --in $db_faa --db $db_prefix";
  
  # 加上 if 判断，如果 .dmnd 文件不存在，才执行 makedb
  if (! -e "$db_prefix.dmnd") {
      print "Building Diamond Database...\n";
      system("$diamond_db");
  }
  
  foreach my $file (@files) {
    my $out = $file;
    $out =~ s/\.$filetype$/.diamond/;  # 避免替换路径中的同名字段
    
    my $evalue  = $config{DIAMOND_EVALUE};
    my $k       = $config{DIAMOND_K};
    my $threads = $config{DIAMOND_THREADS};

    system("$diamond blastx -k $k -e $evalue -p $threads -d $db_prefix -q $file -o $out")
      if $seqtype eq "nucl";
    system("$diamond blastp -k $k -e $evalue -p $threads -d $db_prefix -q $file -o $out")
      if $seqtype eq "prot";
  }
}
elsif ( $method eq "usearch" ) {
  my $usearch = $config{USEARCH_PATH};
  die "Please specify the location of usearch in your config!" if !-e $usearch;
  if ( $filetype =~ /gz/ ) {
    die "Only fastq and fasta files are supported by usearch!\n";
  }
  @files = glob("$workdir/*$filetype");
  foreach my $file (@files) {
    my $out = $file;
    $out =~ s/\.$filetype$/.usearch/;
    
    my $id = $config{USEARCH_ID};
    my $db_faa = $config{DB_FAA};
    system("$usearch -usearch_global $file -db $db_faa -id $id -blast6out $out");
  }
}
elsif ( $method eq "blast" ) {
  my $makeblastdb = $config{MAKEBLASTDB_PATH};
  my $blastp      = $config{BLASTP_PATH};
  my $blastx      = $config{BLASTX_PATH};
  
  die "Please specify the correct paths for BLAST+ in your config!" 
    if !-e $makeblastdb or !-e $blastp or !-e $blastx;

  if ( $filetype =~ /gz|fastq/ ) {
    die "Only fasta files are supported by blast program!";
  }
  @files = glob("$workdir/*$filetype");
  my $db_faa = $config{DB_FAA};
  
  # 如果找不到 BLAST+ 生成的索引文件（.phr, .pin, .psq 等），才执行建库
  if (! -e "$db_faa.pin") {
      system("$makeblastdb -in $db_faa -dbtype prot");
  }
  
  foreach my $file (@files) {
    my $out = $file;
    $out =~ s/\.$filetype$/.blast/;
    
    my $evalue     = $config{BLAST_EVALUE};
    my $max_target = $config{BLAST_MAX_TARGET_SEQS};

    # 使用全新的 BLAST+ 语法
    system("$blastp -db $db_faa -query $file -out $out -outfmt 6 -evalue $evalue -max_target_seqs $max_target")
      if $seqtype eq "prot";
    system("$blastx -db $db_faa -query $file -out $out -outfmt 6 -evalue $evalue -max_target_seqs $max_target")
      if $seqtype eq "nucl";
  }
}

## 5. 读取 Mapping 字典并统计丰度 (Read mapping file and calculate abundance)
my %id2gene;
open( FILE, "$config{MAP_FILE}" ) || die "Cannot open mapping file $config{MAP_FILE}: #1\n";
while (<FILE>) {
  chomp;
  my @items = split( "\t", $_ );
  $id2gene{ $items[0] } = $items[1];
}
close FILE;

my %abundance;
my %samples;
my @sfiles;
@sfiles = glob("$workdir/*diamond") if $method eq "diamond";
@sfiles = glob("$workdir/*usearch") if $method eq "usearch";
@sfiles = glob("$workdir/*blast")   if $method eq "blast";

die "No diamond/usearch/blast output files were detected in $workdir!\n" if $#sfiles==-1;

foreach my $file ( @sfiles ) {
  $file =~ /$workdir\/(.*?)\.$method/;
  my $sample = $1;
  $samples{$sample} = 1;
  my %hit;
  open( FILE, "$file" ) || die "Cannot open output file $file: #2\n";
  while (<FILE>) {
    chomp;
    my @items = split( "\t", $_ );
    my $gene = $id2gene{ $items[1] };
    if ( !$hit{ $items[0] } ) {
      $abundance{$sample}{$gene}++ if $gene;
      $hit{ $items[0] } = 1;
    }
  }
  close FILE;
}

## 6. 读取样本信息及抽平大小 (Read sample info and setup random sampling)
my %size;
my @sizes;
open( FILE, "$sampleinfo" ) || die "Cannot open sample info file $sampleinfo: #3\n";
while (<FILE>) {
  chomp;
  my @items = split( "\t", $_ );
  $size{ $items[0] } = $items[1];
  push( @sizes, $items[1] );
}
close FILE;

foreach my $sample(keys %samples){
  die "$sample was not found in $sampleinfo, please check!\n" if !$size{$sample};
}

@sizes = sort { $a <=> $b } @sizes;

# === 新增逻辑：智能处理抽平 ===
my $rs = $sizes[0]; # 默认取所有样本序列数的最小值
$rs = $randomsampling if defined $randomsampling && $randomsampling ne "";

my %abundance_rs;

# 如果手动设置为 0 或 none，则跳过抽平过程，直接输出绝对丰度
if ( $rs eq "0" || lc($rs) eq "none" ) {
  print "=> Notice: Skipping rarefaction (using absolute counts).\n";
  $rs = "0 (Absolute Counts)"; # 用于在结果文件中做个标记
  # 【修复】：将 $abundance{样本}{基因} 的结构翻转为制表需要的 $abundance_rs{基因}{样本}
  foreach my $s ( keys %abundance ) {
    foreach my $g ( keys %{ $abundance{$s} } ) {
      $abundance_rs{$g}{$s} = $abundance{$s}{$g};
    }
  }
} 
# 【新增修复】：如果用户输入了抽平数值（如 99999 或不填默认取最小），必须调用子函数处理
else {
  print "=> Notice: Random sampling size is set to $rs. Performing rarefaction...\n";
  %abundance_rs = &RandomSampling( \%abundance, \%size, $rs );
}

## 7. 输出最终结果表 (Output final profile)
my @samples = keys %samples;
open( OUT, ">$outfile" ) || die "Cannot write to outfile $outfile: #4\n";
print OUT "#random sampling: $rs\n";
print OUT "Gene\t", join( "\t", @samples ), "\n";
foreach my $gene ( sort keys %abundance_rs ) {
  print OUT "$gene";
  foreach my $sample (@samples) {
    print OUT "\t$abundance_rs{$gene}{$sample}" if $abundance_rs{$gene}{$sample};
    print OUT "\t0"                             if !$abundance_rs{$gene}{$sample};
  }
  print OUT "\n";
}
close OUT;

# 子函数：随机抽样
sub RandomSampling {
  my ( $abundance, $size, $rs ) = @_;
  my %abundance = %$abundance;
  my %size      = %$size;
  my %sum;
  foreach my $sample ( keys %abundance ) {
    foreach my $gene ( keys %{ $abundance{$sample} } ) {
      $sum{$sample} += $abundance{$sample}{$gene};
    }
  }
  my %abundance_rs;
  foreach my $sample ( keys %size ) {
    my $actual_size = $size{$sample}; # 获取当前样本实际拥有的总序列数
    my $target_rs = $rs;              # 获取用户设定的抽平目标值（如 99999）
    # 【核心修复】：防止数组越界崩溃
    # 如果用户设定的抽平值（比如99999）大于样本实际序列数（比如5000），
    # 则强制将抽平上限锁定为实际序列数（5000），防止数组强行切片产生大量 undef 空值
    $target_rs = $actual_size if $target_rs > $actual_size;
    my @array = shuffle( 1 .. $actual_size );
    @array = @array[ 0 .. $target_rs - 1 ]; # 现在最高切到 actual_size，绝对安全
    @array = grep { $_ <= $sum{$sample} } @array;
    my $i = 1;
    foreach my $gene ( keys %{ $abundance{$sample} } ) {
      my @tmp = grep { $_ >= $i && $_ <= ( $abundance{$sample}{$gene} + $i - 1 ) } @array;
      $abundance_rs{$gene}{$sample} = @tmp;
      $i += $abundance{$sample}{$gene};
    }
  }
  return %abundance_rs;
}

sub PrintHelp {
  print "="x50 . "\n";
  print "Usage: perl Ngenes_2606.PL -c <config_file>\n\n";
  print "Note: This is the backend Perl engine for the NCycX pipeline.\n";
  print "Please run the pipeline using the Python wrapper (NCycX_cli.py) \n";
  print "which provides full support for HMM filtering and parameter overrides.\n";
  print "="x50 . "\n";
}
