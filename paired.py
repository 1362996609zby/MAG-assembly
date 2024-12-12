import os
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO # type: ignore
import json
##conda activate metagenomics ##python /mnt/f/binning/paired.py
def run_command(command, log_file):
    """Utility function to run a shell command."""
    with open(log_file, 'a') as log:
        log.write(f"Running command: {command}\n")
        print(f"Running command: {command}")
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            log.write(f"Error running command: {command}\n")
            log.write(result.stderr + '\n')
            print(f"Error running command: {command}")
            print(result.stderr)
        log.write(result.stdout + '\n')
        print(result.stdout)
        return result.stdout

def ensure_empty_dir(directory):
    """Ensure the directory is empty."""
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

def collect_fasta_files(directory, suffix=".fasta"):
    """Collect all fasta files in the given directory."""
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(suffix)]

def convert_concoct_to_fasta(concoct_clustering_file, contigs_fasta_file, output_dir):
    """Convert CONCOCT clustering results to fasta files."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 读取 CONCOCT 聚类结果文件
    clustering_df = pd.read_csv(concoct_clustering_file)
    
    # 读取 contigs 序列
    contigs = SeqIO.to_dict(SeqIO.parse(contigs_fasta_file, "fasta"))
    
    # 创建一个字典来存储每个 bin 的 contigs
    bins = {}
    
    for index, row in clustering_df.iterrows():
        contig_id = row['contig_id']
        bin_id = row['cluster_id']
        if bin_id not in bins:
            bins[bin_id] = []
        bins[bin_id].append(contigs[contig_id])
    
    # 将每个 bin 写入单独的 fasta 文件
    for bin_id, sequences in bins.items():
        bin_fasta_file = os.path.join(output_dir, f"concoct_bin_{bin_id}.fasta")
        SeqIO.write(sequences, bin_fasta_file, "fasta")
    
    return [os.path.join(output_dir, f"concoct_bin_{bin_id}.fasta") for bin_id in bins]

def create_das_tool_input_file(maxbin2_files, metabat2_files, concoct_files, unique_contig_map, output_dir):
    das_tool_input_file = os.path.join(output_dir, "das_tool_input.tsv")
    with open(das_tool_input_file, 'w') as f:
        for bin_files, bin_prefix in [(maxbin2_files, 'maxbin2'), (metabat2_files, 'metabat2'), (concoct_files, 'concoct')]:
            for fasta_file in bin_files:
                bin_id = os.path.basename(fasta_file).replace('.fasta', '').replace('.fa', '')
                with open(fasta_file, 'r') as in_f:
                    for record in SeqIO.parse(in_f, 'fasta'):
                        if (bin_prefix, record.id) in unique_contig_map:
                            unique_contig_id = unique_contig_map[(bin_prefix, record.id)]
                            if unique_contig_id.endswith(f"_{bin_id}"):
                                f.write(f"{unique_contig_id}\t{bin_id}\n")
    return das_tool_input_file

def create_combined_contigs_file(maxbin2_files, metabat2_files, concoct_files, output_file):
    """Create a combined contigs file with unique contig names."""
    # 确保如果文件存在就删除它
    if os.path.exists(output_file):
        os.remove(output_file)

    unique_contig_map = {}
    unique_contig_ids = set()
    with open(output_file, 'w') as out_f:
        for bin_files, prefix in [(maxbin2_files, 'maxbin2'), (metabat2_files, 'metabat2'), (concoct_files, 'concoct')]:
            for fasta_file in bin_files:
                bin_id = os.path.basename(fasta_file).replace('.fasta', '').replace('.fa', '')
                with open(fasta_file, 'r') as in_f:
                    for record in SeqIO.parse(in_f, 'fasta'):
                        unique_id = f"{record.id}_{bin_id}"
                        if unique_id not in unique_contig_ids:
                            unique_contig_ids.add(unique_id)
                            unique_contig_map[(prefix, record.id)] = unique_id
                            record.id = unique_id
                            record.description = ""
                            SeqIO.write(record, out_f, 'fasta')
    return unique_contig_map

def filter_checkm_results(checkm_output_dir, completeness_threshold=80, contamination_threshold=10):
    checkm_results_file = os.path.join(checkm_output_dir, "bin_stats_ext.tsv")

    # 读取文件并解析每一行的字典
    data = []
    with open(checkm_results_file, 'r') as file:
        for line in file:
            bin_name, bin_stats = line.strip().split('\t')
            bin_stats = json.loads(bin_stats.replace("'", "\""))  # 将单引号替换为双引号，并解析为字典
            bin_stats['Bin Name'] = bin_name  # 添加 Bin Name 列
            data.append(bin_stats)

    # 将数据转换为 DataFrame
    checkm_results = pd.DataFrame(data)

    # 过滤结果
    filtered_checkm_results = checkm_results[
        (checkm_results['Completeness'] > completeness_threshold) &
        (checkm_results['Contamination'] < contamination_threshold)
    ]

    # 获取基因组路径
    bins_dir = os.path.join(checkm_output_dir, "bins")
    genome_paths = []

    for bin_name in filtered_checkm_results['Bin Name']:
        # 构建每个 bin 的 genes.fna 文件路径
        original_genome_path = os.path.join(bins_dir, bin_name, "genes.fna")
        renamed_genome_path = os.path.join(bins_dir, bin_name, f"{bin_name}.fna")  # 将文件重命名为 bin_name.fna
        
        if os.path.exists(original_genome_path):
            # 检查目标文件是否已经存在
            if not os.path.exists(renamed_genome_path):
                # 重命名文件
                os.rename(original_genome_path, renamed_genome_path)
            
            # 添加路径到 genome_paths
            genome_paths.append(renamed_genome_path)
        
        elif os.path.exists(renamed_genome_path):
            # 如果已经重命名过，直接使用重命名的路径
            genome_paths.append(renamed_genome_path)
        
        else:
            print(f"Warning: Genome file not found for bin '{bin_name}' at {original_genome_path}")
            genome_paths.append(None)

    # 检查路径数量
    print(f"Number of filtered bins: {len(filtered_checkm_results)}")
    print(f"Number of genome paths: {len(genome_paths)}")

    # 添加基因组路径信息
    filtered_checkm_results['Genome Path'] = genome_paths

    # 保存过滤结果
    filtered_checkm_results_file = os.path.join(checkm_output_dir, "filtered_checkm_results.tsv")
    filtered_checkm_results.to_csv(filtered_checkm_results_file, sep="\t", index=False)

    return filtered_checkm_results_file

def extract_bins_from_das_tool_output(das_tool_output_dir, output_dir):
    """Extract individual bins from DAS Tool output."""
    bin_contig_map_file = os.path.join(das_tool_output_dir, "das_tool_DASTool_contig2bin.tsv")
    if not os.path.exists(bin_contig_map_file):
        raise FileNotFoundError(f"{bin_contig_map_file} not found.")
    
    # 读取文件并指定列名
    bin_contig_map = pd.read_csv(bin_contig_map_file, sep='\t', header=None, names=['contig', 'bin'])
    
    proteins_file = os.path.join(das_tool_output_dir, "das_tool_proteins.faa")
    proteins = SeqIO.to_dict(SeqIO.parse(proteins_file, "fasta"))
    
    # 打印可用的contig ID和从contig2bin文件中提取的ID
    print("Available contig IDs:", list(proteins.keys())[:10])
    print("Contig IDs from das_tool_DASTool_contig2bin.tsv:", bin_contig_map['contig'].unique()[:10])
    
    bins = bin_contig_map.groupby('bin')['contig'].apply(list).to_dict()
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    bin_files = []
    for bin_id, contig_ids in bins.items():
        bin_file = os.path.join(output_dir, f"{bin_id}.faa")
        with open(bin_file, 'w') as f:
            for contig_id in contig_ids:
                # 去掉末尾的数字后缀
                base_contig_id = contig_id.rsplit('_', 1)[0]
                # 检查是否存在于字典中
                matches = [key for key in proteins.keys() if key.startswith(base_contig_id)]
                if matches:
                    for match in matches:
                        SeqIO.write(proteins[match], f, "fasta")
                else:
                    print(f"Warning: {contig_id} not found in proteins")
        bin_files.append(bin_file)
    
    return bin_files

def main():
    log_file = "/mnt/f/binning/pipeline.log"
    #input_file = "/mnt/f/binning/BK2/BK2.scaftigs.fastq"
    #input_files = ["/mnt/f/binning/BK2/BK2.scaftigs.fastq", "/mnt/f/binning/BK2/BK2.scaftigs.fastq"]  # Example list of input files
    #input_files_str = ' '.join(input_files)
    # Paired-end reads paths
    forward_reads = "/mnt/f/binning/RHDCA6/RHDCA6_L1_1.fq"
    reverse_reads = "/mnt/f/binning/RHDCA6/RHDCA6_L1_2.fq"
    
    filtered_fastq_R1 = "/mnt/f/binning/RHDCA6/filtered.fastq_R1.fastq"
    filtered_fastq_R2 = "/mnt/f/binning/RHDCA6/filtered.fastq_R2.fastq"
    
    reference_genome = "/mnt/f/binning/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_sm.toplevel.fa"
    cleaned_fastq_R1 = "/mnt/f/binning/RHDCA6/cleaned_R1.fastq"
    cleaned_fastq_R2 = "/mnt/f/binning/RHDCA6/cleaned_R2.fastq"
    
    megahit_output = "/mnt/f/binning/RHDCA6/megahit_output"
    quast_output = "/mnt/f/binning/RHDCA6/quast_output"
    coverage_file = "/mnt/f/binning/RHDCA6/coverage.tsv"
    binning_output = "/mnt/f/binning/RHDCA6/binning"
    contig_index_prefix = "/mnt/f/binning/RHDCA6/megahit_output/final.contigs"

    ## Create directories if not exist
    os.makedirs(os.path.dirname(filtered_fastq_R1), exist_ok=True)
    os.makedirs(os.path.dirname(cleaned_fastq_R1), exist_ok=True)
    #os.makedirs(megahit_output, exist_ok=True)
    os.makedirs(quast_output, exist_ok=True)
    os.makedirs(binning_output, exist_ok=True)

    # # Step 1: Quality filtering for paired-end reads
    #run_command(f"fastp -i {forward_reads} -I {reverse_reads} -o {filtered_fastq_R1} -O {filtered_fastq_R2}", log_file)

    # # Step 2: Remove host, food, and human sequences
    aligned_sam = "/mnt/f/binning/RHDCA6/aligned.sam"
    aligned_bam = "/mnt/f/binning/RHDCA6/aligned.bam"
    sorted_bam = "/mnt/f/binning/RHDCA6/sorted.bam"
    reference_index = "/mnt/f/binning/RHDCA6/reference_index"

    #run_command(f"bowtie2-build {reference_genome} {reference_index} -p 16", log_file)
    #run_command(f"bowtie2 -x {reference_index} -1 {filtered_fastq_R1} -2 {filtered_fastq_R2} -S {aligned_sam} -p 16", log_file)
    #run_command(f"samtools view -bS {aligned_sam} > {aligned_bam} -@ 16", log_file)
    #run_command(f"samtools sort {aligned_bam} -o {sorted_bam} -@ 16", log_file)
    #run_command(f"samtools index {sorted_bam} -@ 16", log_file)
    #run_command(f"bedtools bamtofastq -i {sorted_bam} -fq {cleaned_fastq_R1} -fq2 {cleaned_fastq_R2}", log_file)


    # # Step 3: Assembly of high-quality reads
    #run_command(f"megahit -1 {cleaned_fastq_R1} -2 {cleaned_fastq_R2} -o {megahit_output} -t 16", log_file)
    #run_command(f"quast {megahit_output}/final.contigs.fa -o {quast_output} --threads 16", log_file)

    # # Step 4: Create Bowtie 2 index for contigs
    #run_command(f"bowtie2-build {megahit_output}/final.contigs.fa {contig_index_prefix}", log_file)

    # # Step 5: Align reads to contigs and calculate coverage
    contig_aligned_sam = "/mnt/f/binning/RHDCA6/contig_aligned.sam"
    contig_aligned_bam = "/mnt/f/binning/RHDCA6/contig_aligned.bam"
    contig_aligned_sorted_bam = "/mnt/f/binning/RHDCA6/contig_aligned_sorted.bam"

    #run_command(f"bowtie2 -x {contig_index_prefix} -1 {cleaned_fastq_R1} -2 {cleaned_fastq_R2} -S {contig_aligned_sam} -p 16", log_file)
    #run_command(f"samtools view -bS {contig_aligned_sam} > {contig_aligned_bam} -@ 16", log_file)
    #run_command(f"samtools sort {contig_aligned_bam} -o {contig_aligned_sorted_bam} -@ 16", log_file)
    #run_command(f"samtools index {contig_aligned_sorted_bam} -@ 16", log_file)
    #run_command(f"samtools idxstats {contig_aligned_sorted_bam} > {coverage_file} -@ 16", log_file)

    
    # 确保目录存在
    os.makedirs(binning_output, exist_ok=True)
    
    # 读取并处理 coverage 文件
    coverage_file = "/mnt/f/binning/RHDCA6/coverage.tsv"
    coverage_file_checked = "/mnt/f/binning/RHDCA6/coverage_checked.tsv"
    #coverage_df = pd.read_csv(coverage_file, sep='\t', header=None, names=['contig', 'length', 'mapped', 'unmapped'])
    #coverage_df['coverage'] = coverage_df['mapped'] / coverage_df['length']
    #coverage_df = coverage_df[['contig', 'coverage']]
    #coverage_df.to_csv(coverage_file_checked, sep='\t', index=False)
    #coverage_df = pd.read_csv(coverage_file_checked, sep='\t')
    
    # 确保列名为字符串
    #coverage_df.columns = coverage_df.columns.astype(str)
    
    # 移除所有非字符串的 contig 行
    #coverage_df = coverage_df[coverage_df['contig'].apply(lambda x: isinstance(x, str))]
    
    # 将所有数值转换为字符串
    #coverage_df = coverage_df.astype(str)
    
    # 检查并移除任何无效行（如空行或包含无效字符的行）
    #coverage_df = coverage_df[coverage_df['contig'].str.strip() != '']
    # Remove rows with missing values
    #coverage_df = coverage_df.dropna()
    
    # 保存处理后的文件
    #coverage_file_checked_corrected = "/mnt/f/binning/RHDCA6/coverage_checked_corrected.tsv"
    #coverage_df.to_csv(coverage_file_checked_corrected, sep='\t', index=False)
    
    # Step 6: Binning
    #run_command(f"run_MaxBin.pl -contig {megahit_output}/final.contigs.fa -out {binning_output}/maxbin2_bin -abund {coverage_file_checked_corrected} -thread 16", log_file)
    #run_command(f"metabat2 -i {megahit_output}/final.contigs.fa -o {binning_output}/metabat2_bins --numThreads 16", log_file)   
    #run_command(f"concoct --composition_file {megahit_output}/final.contigs.fa --coverage_file {coverage_file_checked_corrected} -b {binning_output}/concoct_bins --threads 16", log_file)
 
    # maxbin2, metabat2_bins 和 concoct_bins 的目录
    maxbin2_bin_dir = binning_output
    metabat2_bins_dir = binning_output
    concoct_bins_dir = binning_output
    combined_contigs_file = os.path.join(binning_output, "combined_contigs.fa")

    # 检查路径和文件是否存在
    #if not os.path.exists(maxbin2_bin_dir):
        #print(f"Error: MaxBin2 directory does not exist at path: {maxbin2_bin_dir}")
        #return
    #if not os.path.exists(metabat2_bins_dir):
        #print(f"Error: MetaBat2 bins directory does not exist at path: {metabat2_bins_dir}")
        #return
    #if not os.path.exists(concoct_bins_dir):
       #print(f"Error: CONCOCT bins directory does not exist at path: {concoct_bins_dir}")
        #return

    # 收集 fasta 文件
    #maxbin2_fasta_files = collect_fasta_files(maxbin2_bin_dir, suffix=".fasta")
    #metabat2_fasta_files = collect_fasta_files(metabat2_bins_dir, suffix=".fa")

    # 将 CONCOCT 聚类结果文件转换成 fasta 文件
    #concoct_clustering_file = os.path.join(concoct_bins_dir, "concoct_bins_clustering_gt1000.csv")
    #concoct_fasta_files = convert_concoct_to_fasta(concoct_clustering_file, os.path.join(megahit_output, "final.contigs.fa"), concoct_bins_dir)

    # 创建 combined_contigs_file
    #combined_contigs_file = os.path.join(binning_output, "combined_contigs.fa")
    #unique_contig_map = create_combined_contigs_file(maxbin2_fasta_files, metabat2_fasta_files, concoct_fasta_files, combined_contigs_file)

    # 创建 DAS_Tool 输入文件
    #das_tool_input_file = create_das_tool_input_file(maxbin2_fasta_files, metabat2_fasta_files, concoct_fasta_files, unique_contig_map, binning_output)
    
    #进行这一步之前请手动删除combined_contigs.fa文件
    # Step 7: Integrate MAGs
    #das_tool_output_dir = os.path.join(binning_output, "das_tool")
    #ensure_empty_dir(das_tool_output_dir)  
    #run_command(f"DAS_Tool -i {das_tool_input_file} -c {combined_contigs_file} -o {das_tool_output_dir} --threads 16 --write_bins", log_file)
# 在这之后，手动将dastool生成的11个文件放入das_tool文件夹，使用核酸序列时，这一步不运行
    # Extract individual bins from DAS Tool output
    #extracted_bins_dir = os.path.join(das_tool_output_dir, "extracted_bins")
    #ensure_empty_dir(extracted_bins_dir)
    #bin_files = extract_bins_from_das_tool_output(das_tool_output_dir, extracted_bins_dir)

   # Log the extracted bin files
    #with open(log_file, 'a') as log:
        #log.write(f"Extracted bin files: {bin_files}\n")

#当电脑的运行内存小于40g的时候，使用 --reduced_tree模式
    # Step 8: Quality assessment with CheckM
    extracted_bins_dir = os.path.join(binning_output, "das_tool_DASTool_bins")
    checkm_output_dir = os.path.join(binning_output, "checkm")

    # 运行 CheckM 分析（确保已正确安装并配置 CheckM）
    #run_command(f"checkm lineage_wf --reduced_tree --nt -t 16 -x fa {extracted_bins_dir} {checkm_output_dir}", log_file)

    # 确保 bin_stats_ext.tsv 文件的路径正确
    bin_stats_ext_path = os.path.join(checkm_output_dir, "storage")

    # 筛选符合条件的 MAGs
    if os.path.exists(bin_stats_ext_path):
        filtered_checkm_results_file = filter_checkm_results(bin_stats_ext_path)
    else:
        print(f"File not found: {bin_stats_ext_path}")
        return
#运行这一步之前手动将bins文件夹放入storage文件夹中
    # Step 9: Dereplication with dRep
    genome_list_path = f"{binning_output}/genome_list.txt"

    # 将合格的基因组路径写入 genome_list.txt
    filtered_results = pd.read_csv(filtered_checkm_results_file, sep="\t")
    with open(genome_list_path, 'w') as f:
        for path in filtered_results['Genome Path'].dropna():
            f.write(f"{path}\n")

    run_command(f"dRep dereplicate {binning_output}/drep -g {genome_list_path} --ignoreGenomeQuality", log_file)


    # Step 10: Annotation

    # 设置 GTDB 数据库路径为环境变量
    #os.environ['GTDBTK_DATA_PATH'] = '/mnt/e/zby/binning/GTDB-Tk reference data'
    #run_command(f"gtdbtk classify_wf --genome_dir {binning_output}/drep --out_dir {binning_output}/gtdbtk", log_file)

    # Step 11: K number assignment
    #run_command(f"exec_annotation {megahit_output}/final.contigs.fa -o {binning_output}/kofamscan --cpu 4 -k /path/to/ko_list -p /path/to/profile_hmm", log_file)

if __name__ == "__main__":
    main()


