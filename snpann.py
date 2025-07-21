# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import numpy as np
from collections import defaultdict
import os
from datetime import datetime
from Bio import SeqIO
import re

class VariantAnnotationTool:
    def __init__(self, parent):
        self.root = parent
        
        # 数据存储
        self.vcf_data = None      # 存储VCF变异信息
        self.gwas_data = None     # 存储GWAS结果(P值等)
        self.gff_data = None      # 存储基因注释
        self.fasta_data = None    # 存储参考基因组
        self.annotation_results = None
        
        # 默认配置
        self.config = {
            'p_threshold': 5e-8,
            'output_dir': os.path.expanduser("~")
        }
        
        # 初始化UI
        self.setup_ui()
        
    def setup_ui(self):
        # 主框架
        main_frame = tk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 顶部控制面板
        control_frame = tk.Frame(main_frame)
        control_frame.pack(fill=tk.X, pady=5)
        
        # 数据加载按钮
        btn_load_vcf = tk.Button(control_frame, text="加载VCF文件", command=self.load_vcf)
        btn_load_vcf.pack(side=tk.LEFT, padx=5)
        
        btn_load_gwas = tk.Button(control_frame, text="加载GWAS结果", command=self.load_gwas)
        btn_load_gwas.pack(side=tk.LEFT, padx=5)
        
        btn_load_fasta = tk.Button(control_frame, text="加载参考基因组", command=self.load_fasta)
        btn_load_fasta.pack(side=tk.LEFT, padx=5)
        
        btn_load_gff = tk.Button(control_frame, text="加载基因注释", command=self.load_gff)
        btn_load_gff.pack(side=tk.LEFT, padx=5)
        
        # 阈值设置
        threshold_frame = tk.Frame(control_frame)
        threshold_frame.pack(side=tk.LEFT, padx=10)
        
        tk.Label(threshold_frame, text="P值阈值:").pack(side=tk.LEFT)
        self.p_threshold_var = tk.StringVar(value=str(self.config['p_threshold']))
        tk.Entry(threshold_frame, textvariable=self.p_threshold_var, width=10).pack(side=tk.LEFT)
        
        # 分析按钮
        btn_merge = tk.Button(control_frame, text="合并数据", command=self.merge_data)
        btn_merge.pack(side=tk.LEFT, padx=5)
        
        btn_annotate = tk.Button(control_frame, text="执行注释", command=self.annotate_variants)
        btn_annotate.pack(side=tk.LEFT, padx=5)
        
        # 结果显示区域
        result_frame = tk.Frame(main_frame)
        result_frame.pack(fill=tk.BOTH, expand=True)
        
        # 表格视图
        tree_frame = tk.Frame(result_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        scroll_y = tk.Scrollbar(tree_frame)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        
        scroll_x = tk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
        scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.result_tree = ttk.Treeview(tree_frame, yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        self.result_tree.pack(fill=tk.BOTH, expand=True)
        
        scroll_y.config(command=self.result_tree.yview)
        scroll_x.config(command=self.result_tree.xview)
        
        # 底部按钮
        btn_frame = tk.Frame(main_frame)
        btn_frame.pack(fill=tk.X, pady=5)
        
        btn_export = tk.Button(btn_frame, text="导出结果", command=self.export_results)
        btn_export.pack(side=tk.LEFT, padx=5)
        
        # 状态栏
        self.status_var = tk.StringVar(value="准备就绪")
        status_bar = tk.Label(main_frame, textvariable=self.status_var, bd=1, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X)
    
    def load_vcf(self):
        """加载VCF文件"""
        path = filedialog.askopenfilename(
            initialdir=self.config['output_dir'],
            filetypes=[("VCF文件", "*.vcf *.vcf.gz"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            # 读取VCF文件
            vcf_reader = self.parse_vcf(path)
            self.vcf_data = pd.DataFrame(vcf_reader)
            
            self.status_var.set(f"成功加载VCF文件: {len(self.vcf_data)} 个变异")
            messagebox.showinfo("成功", "VCF文件加载成功")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载VCF文件失败:\n{str(e)}")
            self.status_var.set("加载VCF文件失败")
    
    def parse_vcf(self, vcf_path):
        """解析VCF文件"""
        vcf_data = []
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1])
                var_id = parts[2]
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                filter_ = parts[6]
                info = parts[7]
                
                # 解析INFO字段
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = True
                
                vcf_data.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'ID': var_id,
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': qual,
                    'FILTER': filter_,
                    **info_dict
                })
        
        return vcf_data
    
    def load_gwas(self):
        """加载GWAS结果文件"""
        path = filedialog.askopenfilename(
            initialdir=self.config['output_dir'],
            filetypes=[("CSV文件", "*.csv"), ("文本文件", "*.txt"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            # 尝试自动检测分隔符
            with open(path, 'r') as f:
                first_line = f.readline()
                sep = ',' if ',' in first_line else '\t'
            
            # 读取文件
            self.gwas_data = pd.read_csv(path, sep=sep)
            
            # 自动检测列名
            self.detect_gwas_columns()
            
            self.status_var.set(f"成功加载GWAS结果: {len(self.gwas_data)} 个位点")
            messagebox.showinfo("成功", "GWAS结果加载成功")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载GWAS结果失败:\n{str(e)}")
            self.status_var.set("加载GWAS结果失败")
    
    def detect_gwas_columns(self):
        """自动检测GWAS结果中的必要列"""
        if self.gwas_data is None:
            return
            
        # 可能的列名变体
        snp_cols = ['SNP', 'rsID', 'rsid', 'ID', 'snp']
        chr_cols = ['CHR', 'CHROM', 'chromosome', 'chr']
        pos_cols = ['POS', 'BP', 'position', 'pos']
        p_cols = ['P', 'PVALUE', 'pvalue', 'p_val']
        
        # 查找匹配的列
        self.gwas_snp_col = next((col for col in snp_cols if col in self.gwas_data.columns), None)
        self.gwas_chr_col = next((col for col in chr_cols if col in self.gwas_data.columns), None)
        self.gwas_pos_col = next((col for col in pos_cols if col in self.gwas_data.columns), None)
        self.gwas_p_col = next((col for col in p_cols if col in self.gwas_data.columns), None)
        
        # 验证是否找到所有必要列
        if not all([self.gwas_snp_col, self.gwas_chr_col, self.gwas_pos_col, self.gwas_p_col]):
            messagebox.showwarning("警告", "未能自动识别所有必要列(SNP, CHR, POS, P)")
    
    def load_fasta(self):
        """加载FASTA参考基因组"""
        path = filedialog.askopenfilename(
            initialdir=self.config['output_dir'],
            filetypes=[("FASTA文件", "*.fa *.fasta *.fa.gz *.fasta.gz"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            # 读取FASTA文件
            self.fasta_data = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
            
            self.status_var.set(f"成功加载参考基因组: {len(self.fasta_data)} 条染色体")
            messagebox.showinfo("成功", "FASTA文件加载成功")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载FASTA文件失败:\n{str(e)}")
            self.status_var.set("加载FASTA文件失败")
    
    def load_gff(self):
        """加载GFF3注释文件"""
        path = filedialog.askopenfilename(
            initialdir=self.config['output_dir'],
            filetypes=[("GFF3文件", "*.gff *.gff3"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            # 解析GFF文件
            self.gff_data = self.parse_gff3(path)
            
            self.status_var.set(f"成功加载GFF注释: {len(self.gff_data['genes'])} 个基因")
            messagebox.showinfo("成功", "GFF注释加载成功")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载GFF文件失败:\n{str(e)}")
            self.status_var.set("加载GFF文件失败")
    
    def parse_gff3(self, gff_path):
        """解析GFF3文件，构建基因结构信息"""
        gene_dict = {}
        transcript_dict = {}
        
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                chrom = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                score = parts[5]
                strand = parts[6]
                phase = parts[7]
                attributes = parts[8]
                
                # 解析属性
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                    elif attr:
                        attr_dict[attr] = True
                
                # 存储基因信息
                if feature_type == 'gene':
                    gene_id = attr_dict.get('ID', '')
                    gene_name = attr_dict.get('Name', gene_id)
                    gene_dict[gene_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'name': gene_name,
                        'transcripts': []
                    }
                
                # 存储转录本信息
                elif feature_type == 'mRNA':
                    transcript_id = attr_dict.get('ID', '')
                    parent = attr_dict.get('Parent', '')
                    if parent in gene_dict:
                        transcript_dict[transcript_id] = {
                            'gene': parent,
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'exons': [],
                            'cds': []
                        }
                        gene_dict[parent]['transcripts'].append(transcript_id)
                
                # 存储外显子和CDS信息
                elif feature_type in ['exon', 'CDS']:
                    parent = attr_dict.get('Parent', '')
                    if parent in transcript_dict:
                        feature = {
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'phase': phase
                        }
                        if feature_type == 'exon':
                            transcript_dict[parent]['exons'].append(feature)
                        else:
                            transcript_dict[parent]['cds'].append(feature)
        
        return {
            'genes': gene_dict,
            'transcripts': transcript_dict
        }
    
    def merge_data(self):
        """合并VCF和GWAS结果数据"""
        if self.vcf_data is None:
            messagebox.showwarning("警告", "请先加载VCF文件")
            return
            
        if self.gwas_data is None:
            messagebox.showwarning("警告", "请先加载GWAS结果")
            return
            
        try:
            # 确保VCF数据有ID列
            if 'ID' not in self.vcf_data.columns:
                messagebox.showwarning("警告", "VCF文件缺少ID列，无法合并")
                return
                
            # 确保GWAS数据有必要的列
            if not all([self.gwas_snp_col, self.gwas_p_col]):
                messagebox.showwarning("警告", "GWAS结果缺少必要列")
                return
                
            # 合并数据 - 基于SNP ID
            merged_data = pd.merge(
                self.vcf_data,
                self.gwas_data[[self.gwas_snp_col, self.gwas_p_col]],
                left_on='ID',
                right_on=self.gwas_snp_col,
                how='left'
            )
            
            # 添加P值列名到配置
            self.p_col = self.gwas_p_col
            
            self.vcf_data = merged_data
            self.status_var.set(f"成功合并数据: {len(self.vcf_data)} 个变异")
            messagebox.showinfo("成功", "数据合并完成")
            
        except Exception as e:
            messagebox.showerror("错误", f"数据合并失败:\n{str(e)}")
            self.status_var.set("数据合并失败")
    
    def annotate_variants(self):
        """执行完整的变异注释流程"""
        if self.vcf_data is None:
            messagebox.showwarning("警告", "请先加载并合并VCF和GWAS数据")
            return
            
        if self.fasta_data is None:
            messagebox.showwarning("警告", "请先加载参考基因组FASTA文件")
            return
            
        if self.gff_data is None:
            messagebox.showwarning("警告", "请先加载GFF3注释文件")
            return
            
        try:
            p_threshold = float(self.p_threshold_var.get()) if hasattr(self, 'p_col') and self.p_col else None
            
            # 筛选显著变异
            if p_threshold and self.p_col in self.vcf_data.columns:
                significant = self.vcf_data[self.vcf_data[self.p_col] <= p_threshold].copy()
            else:
                significant = self.vcf_data.copy()
                messagebox.showwarning("注意", "未使用P值筛选，将分析所有变异")
            
            if len(significant) == 0:
                messagebox.showinfo("结果", "未找到符合条件的变异位点")
                return
                
            # 添加注释列
            annotation_cols = [
                'gene', 'transcript', 'feature_type', 
                'amino_acid_change', 'mutation_type', 'protein_effect'
            ]
            for col in annotation_cols:
                if col not in significant.columns:
                    significant[col] = ''
            
            # 进行注释
            for idx, row in significant.iterrows():
                chrom = str(row['CHROM'])
                pos = int(row['POS'])
                ref = row['REF']
                alt = row['ALT']
                
                # 1. 基因水平注释
                gene_info = self.annotate_gene_level(chrom, pos)
                if gene_info:
                    significant.at[idx, 'gene'] = gene_info['gene_name']
                    significant.at[idx, 'feature_type'] = gene_info['feature_type']
                
                # 2. 转录本水平注释
                if gene_info and gene_info['transcript']:
                    transcript_id = gene_info['transcript']
                    transcript_info = self.gff_data['transcripts'].get(transcript_id, {})
                    
                    # 3. 编码区变异分析
                    if transcript_info.get('cds') and gene_info['feature_type'] == 'CDS':
                        protein_effect = self.analyze_protein_effect(
                            chrom, pos, ref, alt, 
                            transcript_info['strand'],
                            transcript_info['cds']
                        )
                        
                        if protein_effect:
                            significant.at[idx, 'transcript'] = transcript_id
                            significant.at[idx, 'amino_acid_change'] = protein_effect['aa_change']
                            significant.at[idx, 'mutation_type'] = protein_effect['mutation_type']
                            significant.at[idx, 'protein_effect'] = protein_effect['effect']
            
            # 显示结果
            self.display_results(significant)
            self.annotation_results = significant
            
            self.status_var.set(f"完成 {len(significant)} 个变异的注释")
            messagebox.showinfo("完成", "变异注释完成")
            
        except Exception as e:
            messagebox.showerror("错误", f"注释过程中出错:\n{str(e)}")
            self.status_var.set("注释失败")
    
    def annotate_gene_level(self, chrom, pos):
        """基因水平注释"""
        if not self.gff_data or 'genes' not in self.gff_data:
            return None
            
        for gene_id, gene in self.gff_data['genes'].items():
            if gene['chrom'] == chrom and gene['start'] <= pos <= gene['end']:
                # 检查是否在转录本区域内
                feature_type = 'gene'
                transcript_id = None
                
                for t_id in gene['transcripts']:
                    transcript = self.gff_data['transcripts'].get(t_id, {})
                    if transcript.get('start') <= pos <= transcript.get('end'):
                        feature_type = 'transcript'
                        transcript_id = t_id
                        
                        # 检查是否在外显子或CDS区域
                        for exon in transcript.get('exons', []):
                            if exon['start'] <= pos <= exon['end']:
                                feature_type = 'exon'
                                break
                                
                        for cds in transcript.get('cds', []):
                            if cds['start'] <= pos <= cds['end']:
                                feature_type = 'CDS'
                                break
                                
                        break
                
                return {
                    'gene_id': gene_id,
                    'gene_name': gene['name'],
                    'feature_type': feature_type,
                    'transcript': transcript_id,
                    'strand': gene['strand']
                }
        
        return None
    
    def analyze_protein_effect(self, chrom, pos, ref, alt, strand, cds_regions):
        """分析蛋白质水平效应"""
        try:
            # 获取参考序列
            if chrom not in self.fasta_data:
                return None
                
            ref_seq = str(self.fasta_data[chrom].seq)
            
            # 确定变异在CDS中的位置
            cds_pos, codon_num, codon_pos = self.locate_in_cds(pos, cds_regions, strand)
            if cds_pos is None:
                return None
                
            # 获取参考密码子
            ref_codon = self.get_reference_codon(ref_seq, pos, cds_pos, codon_num, codon_pos, strand)
            if not ref_codon:
                return None
                
            # 获取变异后密码子
            alt_codon = self.get_alternate_codon(ref_codon, ref, alt, codon_pos, strand)
            
            # 翻译密码子
            ref_aa = self.translate_codon(ref_codon)
            alt_aa = self.translate_codon(alt_codon)
            
            # 确定突变类型
            mutation_type = self.determine_mutation_type(ref_aa, alt_aa)
            
            # 构建结果
            aa_change = f"{ref_aa}{codon_num+1}{alt_aa}" if ref_aa != alt_aa else "同义突变"
            
            return {
                'aa_change': aa_change,
                'mutation_type': mutation_type,
                'effect': self.get_effect_description(mutation_type, ref_aa, alt_aa)
            }
            
        except Exception as e:
            print(f"蛋白质效应分析错误: {str(e)}")
            return None
    
    def locate_in_cds(self, pos, cds_regions, strand):
        """确定变异在CDS中的位置"""
        # 按基因组坐标排序CDS区域
        if strand == '+':
            sorted_cds = sorted(cds_regions, key=lambda x: x['start'])
        else:
            sorted_cds = sorted(cds_regions, key=lambda x: x['start'], reverse=True)
        
        # 计算CDS中的相对位置
        cds_pos = 0
        for cds in sorted_cds:
            if cds['start'] <= pos <= cds['end']:
                # 在CDS区域内的位置
                offset = pos - cds['start'] if strand == '+' else cds['end'] - pos
                cds_pos += offset
                
                # 计算密码子位置
                codon_num = cds_pos // 3
                codon_pos = cds_pos % 3
                return cds_pos, codon_num, codon_pos
                
            cds_pos += cds['end'] - cds['start'] + 1
        
        return None, None, None
    
    def get_reference_codon(self, ref_seq, pos, cds_pos, codon_num, codon_pos, strand):
        """获取参考密码子"""
        # 确定密码子起始位置
        if strand == '+':
            codon_start = pos - codon_pos
        else:
            codon_start = pos + codon_pos
        
        # 检查边界
        if codon_start < 1 or codon_start + 2 > len(ref_seq):
            return None
            
        # 获取密码子序列
        codon = ref_seq[codon_start-1:codon_start+2].upper()
        
        # 反向互补链需要反转
        if strand == '-':
            codon = self.reverse_complement(codon)
        
        return codon
    
    def get_alternate_codon(self, ref_codon, ref, alt, codon_pos, strand):
        """获取变异后密码子"""
        # 处理SNP
        if len(ref) == 1 and len(alt) == 1:
            alt_codon = list(ref_codon)
            if strand == '+':
                alt_codon[codon_pos] = alt.upper()
            else:
                alt_codon[2-codon_pos] = self.reverse_complement(alt.upper())
            return ''.join(alt_codon)
        
        # 处理插入/缺失 (简化处理)
        elif len(ref) > len(alt):  # 缺失
            return '---'  # 标记为缺失
        elif len(alt) > len(ref):  # 插入
            return '+++'  # 标记为插入
        
        return ref_codon
    
    def translate_codon(self, codon):
        """翻译密码子为氨基酸"""
        codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
        }
        
        if len(codon) != 3:
            return '?'
        
        return codon_table.get(codon, '?')
    
    def reverse_complement(self, seq):
        """获取反向互补序列"""
        comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        return ''.join([comp.get(base, 'N') for base in seq[::-1]])
    
    def determine_mutation_type(self, ref_aa, alt_aa):
        """确定突变类型"""
        if ref_aa == alt_aa:
            return "同义突变"
        elif alt_aa == '*':
            return "无义突变"
        elif ref_aa == '*' and alt_aa != '*':
            return "终止密码子获得"
        elif ref_aa != '*' and alt_aa == '*':
            return "提前终止"
        else:
            return "错义突变"
    
    def get_effect_description(self, mutation_type, ref_aa, alt_aa):
        """获取效应描述"""
        descriptions = {
            "同义突变": "不改变氨基酸序列",
            "错义突变": f"氨基酸改变: {ref_aa}→{alt_aa}",
            "无义突变": "引入终止密码子",
            "提前终止": "导致蛋白质截短",
            "终止密码子获得": "消除终止密码子"
        }
        return descriptions.get(mutation_type, "未知效应")
    
    def display_results(self, df):
        """在表格中显示结果"""
        self.result_tree.delete(*self.result_tree.get_children())
        
        # 设置列
        columns = list(df.columns)
        self.result_tree["columns"] = columns
        self.result_tree["show"] = "headings"
        
        # 设置列宽和标题
        for col in columns:
            self.result_tree.heading(col, text=col)
            self.result_tree.column(col, width=120, anchor=tk.CENTER)
        
        # 添加数据
        for _, row in df.iterrows():
            self.result_tree.insert("", tk.END, values=list(row))
    
    def export_results(self):
        """导出结果"""
        if self.annotation_results is None or len(self.annotation_results) == 0:
            messagebox.showwarning("警告", "没有可导出的结果")
            return
            
        default_name = f"annotated_variants_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        path = filedialog.asksaveasfilename(
            initialdir=self.config['output_dir'],
            initialfile=default_name,
            defaultextension=".csv",
            filetypes=[("CSV文件", "*.csv"), ("Excel文件", "*.xlsx"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            if path.endswith('.csv'):
                self.annotation_results.to_csv(path, index=False)
            else:
                self.annotation_results.to_excel(path, index=False)
                
            self.status_var.set(f"结果已保存到: {path}")
            messagebox.showinfo("成功", f"结果已保存到:\n{path}")
            
        except Exception as e:
            messagebox.showerror("错误", f"保存失败:\n{str(e)}")
            self.status_var.set("保存失败")



if __name__ == "__main__":
    root = tk.Tk()
    app = VariantAnnotationTool(root)
    root.mainloop()

def create_frame(parent):
    frame = tk.Frame(parent)
    VariantAnnotationTool(frame)  
    return frame