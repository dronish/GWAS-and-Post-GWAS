# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, StringVar, BooleanVar, IntVar
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import os
from datetime import datetime
from collections import defaultdict

class MolecularDesignTool:
    def __init__(self, parent):
        self.root = parent
        
        # 数据存储
        self.input_data = None
        self.fasta_data = None
        self.filtered_variants = None
        self.design_results = None
        self.current_design_type = None
        
        # KSAP探针序列
        self.ksap_probes = {
            'FAM': 'GAAGGTGACCAAGTTCATGCT',
            'VIC': 'GAAGGTCGGAGTCAACGGATT'
        }
        
        # 默认配置
        self.config = {
            'p_threshold': 5e-8,
            'output_dir': os.path.expanduser("~"),
            'probe_length': 20,
            'primer_tm': 60.0,
            'product_size_min': 100,
            'product_size_max': 300,
            'primer_size_range': (18, 25),
            'tail_length': 20,
            'probe_tm': 65.0,
            'tm_diff_threshold': 5.0,
            'primer_pairs': 3
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
        btn_load_data = tk.Button(control_frame, text="加载输入文件", command=self.load_input_data)
        btn_load_data.pack(side=tk.LEFT, padx=5)
        
        btn_load_fasta = tk.Button(control_frame, text="加载参考基因组", command=self.load_fasta)
        btn_load_fasta.pack(side=tk.LEFT, padx=5)
        
        # 筛选设置
        filter_frame = tk.Frame(control_frame)
        filter_frame.pack(side=tk.LEFT, padx=10)
        
        tk.Label(filter_frame, text="P值阈值:").pack(side=tk.LEFT)
        self.p_threshold_var = StringVar(value=str(self.config['p_threshold']))
        tk.Entry(filter_frame, textvariable=self.p_threshold_var, width=10).pack(side=tk.LEFT)
        
        # 突变类型筛选
        self.mutation_types = ['同义突变', '错义突变', '提前终止']
        self.mutation_vars = {m: BooleanVar(value=True) for m in self.mutation_types}
        
        mutation_menu = tk.Menubutton(filter_frame, text="突变类型", relief=tk.RAISED)
        mutation_menu.pack(side=tk.LEFT, padx=5)
        mutation_menu.menu = tk.Menu(mutation_menu, tearoff=0)
        mutation_menu["menu"] = mutation_menu.menu
        
        for m in self.mutation_types:
            mutation_menu.menu.add_checkbutton(label=m, variable=self.mutation_vars[m])
        
        # 设计类型选择
        self.design_type_var = StringVar(value="probe")
        design_type_frame = tk.Frame(control_frame)
        design_type_frame.pack(side=tk.LEFT, padx=10)
        
        tk.Radiobutton(design_type_frame, text="探针设计", variable=self.design_type_var, 
                      value="probe", command=self.setup_design_params).pack(side=tk.LEFT)
        tk.Radiobutton(design_type_frame, text="扩增引物", variable=self.design_type_var, 
                      value="primer", command=self.setup_design_params).pack(side=tk.LEFT)
        tk.Radiobutton(design_type_frame, text="KSAP引物", variable=self.design_type_var, 
                      value="kasp", command=self.setup_design_params).pack(side=tk.LEFT)
        
        # 设计参数设置
        self.design_params_frame = tk.Frame(control_frame)
        self.design_params_frame.pack(side=tk.LEFT, padx=10)
        self.setup_design_params()
        
        # 操作按钮
        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(side=tk.LEFT, padx=10)
        
        btn_filter = tk.Button(btn_frame, text="筛选位点", command=self.filter_variants)
        btn_filter.pack(side=tk.LEFT, padx=5)
        
        btn_design = tk.Button(btn_frame, text="执行设计", command=self.execute_design)
        btn_design.pack(side=tk.LEFT, padx=5)
        
        btn_export = tk.Button(btn_frame, text="导出结果", command=self.export_results)
        btn_export.pack(side=tk.LEFT, padx=5)
        
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
        
        self.result_tree = ttk.Treeview(tree_frame, yscrollcommand=scroll_y.set, 
                                      xscrollcommand=scroll_x.set)
        self.result_tree.pack(fill=tk.BOTH, expand=True)
        
        scroll_y.config(command=self.result_tree.yview)
        scroll_x.config(command=self.result_tree.xview)
        
        # 状态栏
        self.status_var = StringVar(value="准备就绪")
        status_bar = tk.Label(main_frame, textvariable=self.status_var, bd=1, 
                            relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X)
    
    def setup_design_params(self):
        """根据设计类型设置参数控件"""
        # 清除现有控件
        for widget in self.design_params_frame.winfo_children():
            widget.destroy()
        
        design_type = self.design_type_var.get()
        
        if design_type == "probe":
            # 探针设计参数
            tk.Label(self.design_params_frame, text="探针长度:").pack(side=tk.LEFT)
            self.probe_len_var = StringVar(value=str(self.config['probe_length']))
            tk.Entry(self.design_params_frame, textvariable=self.probe_len_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="Tm目标:").pack(side=tk.LEFT)
            self.probe_tm_var = StringVar(value=str(self.config['probe_tm']))
            tk.Entry(self.design_params_frame, textvariable=self.probe_tm_var, 
                    width=5).pack(side=tk.LEFT)
            
        elif design_type == "primer":
            # 扩增引物参数
            tk.Label(self.design_params_frame, text="Tm目标:").pack(side=tk.LEFT)
            self.primer_tm_var = StringVar(value=str(self.config['primer_tm']))
            tk.Entry(self.design_params_frame, textvariable=self.primer_tm_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="产物长度范围:").pack(side=tk.LEFT)
            self.product_size_min_var = StringVar(value=str(self.config['product_size_min']))
            tk.Entry(self.design_params_frame, textvariable=self.product_size_min_var, 
                    width=5).pack(side=tk.LEFT)
            tk.Label(self.design_params_frame, text="-").pack(side=tk.LEFT)
            self.product_size_max_var = StringVar(value=str(self.config['product_size_max']))
            tk.Entry(self.design_params_frame, textvariable=self.product_size_max_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="Tm差异阈值:").pack(side=tk.LEFT)
            self.tm_diff_var = StringVar(value=str(self.config['tm_diff_threshold']))
            tk.Entry(self.design_params_frame, textvariable=self.tm_diff_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="引物对数:").pack(side=tk.LEFT)
            self.primer_pairs_var = IntVar(value=self.config['primer_pairs'])
            tk.Spinbox(self.design_params_frame, from_=1, to=5, 
                      textvariable=self.primer_pairs_var, width=3).pack(side=tk.LEFT)
            
        elif design_type == "kasp":
            # KSAP引物参数
            tk.Label(self.design_params_frame, text="Tm目标:").pack(side=tk.LEFT)
            self.primer_tm_var = StringVar(value=str(self.config['primer_tm']))
            tk.Entry(self.design_params_frame, textvariable=self.primer_tm_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="产物长度范围:").pack(side=tk.LEFT)
            self.product_size_min_var = StringVar(value=str(self.config['product_size_min']))
            tk.Entry(self.design_params_frame, textvariable=self.product_size_min_var, 
                    width=5).pack(side=tk.LEFT)
            tk.Label(self.design_params_frame, text="-").pack(side=tk.LEFT)
            self.product_size_max_var = StringVar(value=str(self.config['product_size_max']))
            tk.Entry(self.design_params_frame, textvariable=self.product_size_max_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="尾长:").pack(side=tk.LEFT)
            self.tail_length_var = StringVar(value=str(self.config['tail_length']))
            tk.Entry(self.design_params_frame, textvariable=self.tail_length_var, 
                    width=5).pack(side=tk.LEFT)
            
            tk.Label(self.design_params_frame, text="Tm差异阈值:").pack(side=tk.LEFT)
            self.tm_diff_var = StringVar(value=str(self.config['tm_diff_threshold']))
            tk.Entry(self.design_params_frame, textvariable=self.tm_diff_var, 
                    width=5).pack(side=tk.LEFT)
    
    def load_input_data(self):
        """加载输入文件(GWAS结果或注释结果)"""
        path = filedialog.askopenfilename(
            initialdir=self.config['output_dir'],
            filetypes=[("CSV文件", "*.csv"), ("文本文件", "*.txt"), ("Excel文件", "*.xlsx"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            # 尝试自动检测分隔符
            if path.endswith('.csv'):
                self.input_data = pd.read_csv(path)
            elif path.endswith('.xlsx') or path.endswith('.xls'):
                self.input_data = pd.read_excel(path)
            else:
                with open(path, 'r') as f:
                    first_line = f.readline()
                    sep = ',' if ',' in first_line else '\t'
                self.input_data = pd.read_csv(path, sep=sep)
            
            # 检查必要列
            self.check_required_columns()
            
            self.status_var.set(f"成功加载输入文件: {len(self.input_data)} 行")
            messagebox.showinfo("成功", "输入文件加载成功")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载输入文件失败:\n{str(e)}")
            self.status_var.set("加载输入文件失败")
    
    def check_required_columns(self):
        """检查必要列并提示用户"""
        if self.input_data is None:
            return False
            
        # 检测必要列
        required_cols = {
            'chr': ['CHR', 'CHROM', 'chromosome', 'chr'],
            'pos': ['POS', 'BP', 'position', 'pos'],
            'ref': ['REF', 'ref', 'reference'],
            'alt': ['ALT', 'alt', 'alternate']
        }
        
        # 查找匹配的列
        self.chr_col = self.find_matching_column(required_cols['chr'])
        self.pos_col = self.find_matching_column(required_cols['pos'])
        self.ref_col = self.find_matching_column(required_cols['ref'])
        self.alt_col = self.find_matching_column(required_cols['alt'])
        
        # 检查必要列是否存在
        missing_cols = []
        if not self.chr_col: missing_cols.append("CHR")
        if not self.pos_col: missing_cols.append("POS")
        if not self.ref_col: missing_cols.append("REF")
        if not self.alt_col: missing_cols.append("ALT")
        
        if missing_cols:
            messagebox.showerror("错误", f"缺少必要列: {', '.join(missing_cols)}")
            return False
        
        # 检查可选列并提示用户
        self.p_col = self.find_matching_column(['P', 'PVALUE', 'pvalue', 'p_val'])
        self.mutation_col = self.find_matching_column(['mutation_type', '突变类型', 'effect', '蛋白效应'])
        
        if not self.p_col:
            messagebox.showwarning("提示", "未检测到P值列，将不使用P值筛选")
        
        if not self.mutation_col:
            messagebox.showwarning("提示", "未检测到突变类型列，将不使用突变类型筛选")
        
        return True
    
    def find_matching_column(self, possible_names):
        """在数据框中查找匹配的列名"""
        for name in possible_names:
            if name in self.input_data.columns:
                return name
        return None
    
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
    
    def filter_variants(self):
        """筛选变异位点"""
        if self.input_data is None:
            messagebox.showwarning("警告", "请先加载输入文件")
            return
            
        try:
            # 获取筛选参数
            p_threshold = float(self.p_threshold_var.get()) if self.p_col else None
            selected_mutations = [m for m, var in self.mutation_vars.items() if var.get()] if self.mutation_col else None
            
            # 应用筛选条件
            filtered = self.input_data.copy()
            
            # P值筛选
            if p_threshold and self.p_col:
                filtered = filtered[filtered[self.p_col] <= p_threshold]
            
            # 突变类型筛选
            if selected_mutations and self.mutation_col:
                filtered = filtered[filtered[self.mutation_col].isin(selected_mutations)]
            
            if len(filtered) == 0:
                messagebox.showinfo("结果", "没有符合条件的变异位点")
                return
                
            self.filtered_variants = filtered
            self.display_results(filtered)
            
            self.status_var.set(f"筛选到 {len(filtered)} 个变异位点")
            messagebox.showinfo("成功", f"筛选到 {len(filtered)} 个符合条件的变异位点")
            
        except Exception as e:
            messagebox.showerror("错误", f"筛选过程中出错:\n{str(e)}")
            self.status_var.set("筛选失败")
    
    def execute_design(self):
        """根据选择的设计类型执行设计"""
        if self.filtered_variants is None or len(self.filtered_variants) == 0:
            messagebox.showwarning("警告", "请先筛选变异位点")
            return
            
        if self.fasta_data is None:
            messagebox.showwarning("警告", "请先加载参考基因组FASTA文件")
            return
            
        design_type = self.design_type_var.get()
        self.current_design_type = design_type
        
        try:
            if design_type == "probe":
                self.design_probes()
            elif design_type == "primer":
                self.design_primers()
            elif design_type == "kasp":
                self.design_kasp_primers()
            
            messagebox.showinfo("成功", f"{self.get_design_type_name()}设计完成")
            
        except Exception as e:
            messagebox.showerror("错误", f"设计过程中出错:\n{str(e)}")
            self.status_var.set(f"{self.get_design_type_name()}设计失败")
    
    def get_design_type_name(self):
        """获取当前设计类型的名称"""
        names = {
            "probe": "探针",
            "primer": "扩增引物",
            "kasp": "KSAP引物"
        }
        return names.get(self.current_design_type, "设计")
    
    def design_probes(self):
        """设计探针"""
        probe_length = int(self.probe_len_var.get())
        target_tm = float(self.probe_tm_var.get())
        
        # 准备结果存储
        results = []
        
        # 对每个变异位点进行设计
        for _, row in self.filtered_variants.iterrows():
            chrom = str(row[self.chr_col])
            pos = int(row[self.pos_col])
            ref = row[self.ref_col]
            alt = row[self.alt_col]
            
            # 获取参考序列
            if chrom not in self.fasta_data:
                continue
                
            seq_record = self.fasta_data[chrom]
            seq = str(seq_record.seq).upper()
            
            # 设计探针
            probe_info = self.design_probe_sequence(seq, pos, ref, alt, probe_length, target_tm)
            
            # 保存结果
            result = {
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'Probe_Sequence': probe_info['sequence'],
                'Probe_Length': probe_length,
                'Probe_Tm': probe_info['tm'],
                'Probe_Start': probe_info['start'],
                'Probe_End': probe_info['end'],
                'Mutation_Position': probe_info['mut_pos']
            }
            
            if self.p_col:
                result['P_value'] = row[self.p_col]
            if self.mutation_col:
                result['Mutation_Type'] = row[self.mutation_col]
            
            results.append(result)
        
        # 显示结果
        self.design_results = pd.DataFrame(results)
        self.display_results(self.design_results)
        
        self.status_var.set(f"完成 {len(results)} 个位点的探针设计")
    
    def design_probe_sequence(self, seq, pos, ref, alt, length, target_tm):
        """设计探针序列，并突出显示突变位点"""
        # 转换为0-based坐标
        pos = pos - 1
        
        # 确保突变位点在序列范围内
        if pos < 0 or pos >= len(seq):
            return {
                'sequence': '',
                'tm': 0,
                'start': 0,
                'end': 0,
                'mut_pos': 0
            }
        
        # 计算探针起始和结束位置
        half_len = length // 2
        start = max(0, pos - half_len)
        end = min(len(seq), pos + half_len + (length % 2))
        
        # 获取探针序列
        probe_seq = seq[start:end]
        
        # 标记变异位置
        mut_pos_in_probe = pos - start
        if 0 <= mut_pos_in_probe < len(probe_seq):
            # 确保突变位点正确
            if probe_seq[mut_pos_in_probe] == ref:
                # 构建标记序列
                marked_seq = (probe_seq[:mut_pos_in_probe] + 
                            f"[{ref}/{alt}]" + 
                            probe_seq[mut_pos_in_probe+1:])
            else:
                # 如果参考序列不匹配，只显示参考序列
                marked_seq = probe_seq
        else:
            marked_seq = probe_seq
        
        # 计算Tm值
        tm = mt.Tm_NN(Seq(probe_seq))
        
        return {
            'sequence': marked_seq,
            'tm': round(tm, 1),
            'start': start + 1,  # 转换为1-based
            'end': end,
            'mut_pos': pos + 1  # 转换为1-based
        }
    
    def design_primers(self):
        """设计扩增引物"""
        target_tm = float(self.primer_tm_var.get())
        min_size = int(self.product_size_min_var.get())
        max_size = int(self.product_size_max_var.get())
        tm_diff_threshold = float(self.tm_diff_var.get())
        primer_pairs = int(self.primer_pairs_var.get())
        
        # 准备结果存储
        results = []
        
        # 对每个变异位点进行设计
        for _, row in self.filtered_variants.iterrows():
            chrom = str(row[self.chr_col])
            pos = int(row[self.pos_col])
            ref = row[self.ref_col]
            alt = row[self.alt_col]
            
            # 获取参考序列
            if chrom not in self.fasta_data:
                continue
                
            seq_record = self.fasta_data[chrom]
            seq = str(seq_record.seq).upper()
            
            # 设计引物
            primers_list = self.design_multiple_pcr_primers(
                seq, pos, target_tm, min_size, max_size, 
                tm_diff_threshold, primer_pairs)
            
            # 为每对引物创建结果条目
            for i, primers in enumerate(primers_list, 1):
                result = {
                    'CHROM': chrom,
                    'POS': pos,
                    'REF': ref,
                    'ALT': alt,
                    'Primer_Pair': i,
                    'Forward_Primer': primers['F'],
                    'Reverse_Primer': primers['R'],
                    'Product_Size': primers['size'],
                    'F_Tm': primers['F_tm'],
                    'R_Tm': primers['R_tm'],
                    'F_Start': primers['F_start'],
                    'F_End': primers['F_end'],
                    'R_Start': primers['R_start'],
                    'R_End': primers['R_end'],
                    'Variant_Position': primers['var_pos']
                }
                
                if self.p_col:
                    result['P_value'] = row[self.p_col]
                if self.mutation_col:
                    result['Mutation_Type'] = row[self.mutation_col]
                
                results.append(result)
        
        # 显示结果
        self.design_results = pd.DataFrame(results)
        self.display_results(self.design_results)
        
        self.status_var.set(f"完成 {len(results)} 个位点的扩增引物设计 ({primer_pairs}对/位点)")
    
    def design_multiple_pcr_primers(self, seq, pos, target_tm, min_size, max_size, 
                                 tm_diff_threshold, num_pairs):
        """设计多对PCR引物"""
        variant_pos = pos - 1  # 转换为0-based坐标
        primers_found = []
        
        min_len, max_len = self.config['primer_size_range']
        
        # 使用字典记录已找到的引物对，避免重复
        found_combinations = defaultdict(bool)
        
        # 尝试不同的产物大小
        for product_size in range(min_size, max_size + 1):
            # 计算目标区域
            upstream_len = product_size // 2
            downstream_len = product_size - upstream_len
            
            # 获取目标区域序列
            target_start = max(0, variant_pos - upstream_len)
            target_end = min(len(seq), variant_pos + downstream_len)
            
            # 跳过太小或太大的区域
            if target_end - target_start < min_size:
                continue
                
            target_seq = seq[target_start:target_end]
            
            # 尝试不同的引物长度
            for f_len in range(min_len, max_len + 1):
                for r_len in range(min_len, max_len + 1):
                    # 前引物
                    f_primer = target_seq[:f_len]
                    f_tm = mt.Tm_NN(Seq(f_primer))
                    
                    # 后引物
                    r_seq = str(Seq(target_seq[-r_len:]).reverse_complement())
                    r_tm = mt.Tm_NN(Seq(r_seq))
                    
                    # 计算与目标Tm的差异
                    f_diff = abs(f_tm - target_tm)
                    r_diff = abs(r_tm - target_tm)
                    
                    # 检查产物大小
                    product_len = len(target_seq) - f_len - r_len
                    
                    # 检查条件
                    if (f_diff <= tm_diff_threshold and 
                        r_diff <= tm_diff_threshold and 
                        min_size <= product_len <= max_size):
                        
                        # 创建唯一键避免重复
                        primer_key = (f_primer, r_seq)
                        
                        if not found_combinations[primer_key]:
                            found_combinations[primer_key] = True
                            
                            primers_found.append({
                                'F': f_primer,
                                'R': r_seq,
                                'F_tm': round(f_tm, 1),
                                'R_tm': round(r_tm, 1),
                                'size': product_len,
                                'F_start': target_start + 1,  # 1-based
                                'F_end': target_start + f_len,
                                'R_start': target_end - r_len + 1,
                                'R_end': target_end,
                                'var_pos': variant_pos + 1 - target_start  # 在产物中的位置
                            })
                            
                            # 如果已经找到足够的引物对，返回结果
                            if len(primers_found) >= num_pairs:
                                return primers_found
        
        return primers_found[:num_pairs]  # 返回找到的引物对，不超过请求的数量
    
    def design_kasp_primers(self):
        """设计KSAP引物"""
        target_tm = float(self.primer_tm_var.get())
        min_size = int(self.product_size_min_var.get())
        max_size = int(self.product_size_max_var.get())
        tail_length = int(self.tail_length_var.get())
        tm_diff_threshold = float(self.tm_diff_var.get())
        
        # 准备结果存储
        results = []
        
        # 对每个变异位点进行设计
        for _, row in self.filtered_variants.iterrows():
            chrom = str(row[self.chr_col])
            pos = int(row[self.pos_col])
            ref = row[self.ref_col]
            alt = row[self.alt_col]
            
            # 获取参考序列
            if chrom not in self.fasta_data:
                continue
                
            seq_record = self.fasta_data[chrom]
            seq = str(seq_record.seq).upper()
            
            # 设计KSAP引物
            kasp_primers = self.design_kasp_primers_for_variant(
                seq, pos, ref, alt, target_tm, min_size, max_size, tail_length, tm_diff_threshold)
            
            # 保存结果
            result = {
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'FAM_Tail_Primer': kasp_primers['FAM']['tail_primer'],
                'FAM_Allele_Primer': kasp_primers['FAM']['allele_primer'],
                'FAM_Sequence': kasp_primers['FAM']['sequence'],
                'VIC_Tail_Primer': kasp_primers['VIC']['tail_primer'],
                'VIC_Allele_Primer': kasp_primers['VIC']['allele_primer'],
                'VIC_Sequence': kasp_primers['VIC']['sequence'],
                'Common_Reverse_Primer': kasp_primers['common_reverse'],
                'Product_Size': kasp_primers['product_size'],
                'FAM_Tm': kasp_primers['FAM']['tm'],
                'VIC_Tm': kasp_primers['VIC']['tm'],
                'Reverse_Tm': kasp_primers['reverse_tm'],
                'FAM_Tail_Length': len(kasp_primers['FAM']['tail_primer']) - len(self.ksap_probes['FAM']),
                'VIC_Tail_Length': len(kasp_primers['VIC']['tail_primer']) - len(self.ksap_probes['VIC']),
                'Reverse_Length': len(kasp_primers['common_reverse'])
            }
            
            if self.p_col:
                result['P_value'] = row[self.p_col]
            if self.mutation_col:
                result['Mutation_Type'] = row[self.mutation_col]
            
            results.append(result)
        
        # 显示结果
        self.design_results = pd.DataFrame(results)
        self.display_results(self.design_results)
        
        self.status_var.set(f"完成 {len(results)} 个位点的KSAP引物设计")
    
    def design_kasp_primers_for_variant(self, seq, pos, ref, alt, target_tm, min_size, max_size, tail_length, tm_diff_threshold):
        """为单个变异设计KSAP引物"""
        # 计算目标区域
        variant_pos = pos - 1  # 转换为0-based坐标
        
        # 寻找最佳引物组合
        best_primers = None
        best_diff = float('inf')
        
        # 尝试不同的产物大小
        for product_size in range(min_size, max_size + 1):
            # 计算目标区域
            upstream_len = product_size // 2
            downstream_len = product_size - upstream_len
            
            # 获取目标区域序列
            target_start = max(0, variant_pos - upstream_len)
            target_end = min(len(seq), variant_pos + downstream_len)
            
            # 跳过太小或太大的区域
            if target_end - target_start < min_size:
                continue
                
            target_seq = seq[target_start:target_end]
            
            # 调整变异位置坐标
            variant_pos_in_target = variant_pos - target_start
            
            # 设计等位基因特异性引物
            fam_primer = self.design_allele_specific_primer(
                target_seq, variant_pos_in_target, ref, 'FAM', tail_length, target_tm)
            
            vic_primer = self.design_allele_specific_primer(
                target_seq, variant_pos_in_target, alt, 'VIC', tail_length, target_tm)
            
            # 设计通用反向引物
            reverse_primer = self.design_reverse_primer(target_seq, target_tm, tm_diff_threshold)
            
            # 检查Tm差异
            fam_diff = abs(fam_primer['tm'] - target_tm)
            vic_diff = abs(vic_primer['tm'] - target_tm)
            rev_diff = abs(reverse_primer['tm'] - target_tm)
            
            # 计算产物大小
            product_len = len(target_seq) - len(fam_primer['allele_primer']) - len(reverse_primer['sequence'])
            
            # 检查条件
            if (fam_diff <= tm_diff_threshold and 
                vic_diff <= tm_diff_threshold and 
                rev_diff <= tm_diff_threshold and 
                min_size <= product_len <= max_size):
                
                total_diff = fam_diff + vic_diff + rev_diff
                
                if total_diff < best_diff:
                    best_diff = total_diff
                    best_primers = {
                        'FAM': fam_primer,
                        'VIC': vic_primer,
                        'common_reverse': reverse_primer['sequence'],
                        'product_size': product_len,
                        'reverse_tm': reverse_primer['tm']
                    }
        
        return best_primers if best_primers else {
            'FAM': {'tail_primer': '', 'allele_primer': '', 'sequence': '', 'tm': 0},
            'VIC': {'tail_primer': '', 'allele_primer': '', 'sequence': '', 'tm': 0},
            'common_reverse': '',
            'product_size': 0,
            'reverse_tm': 0
        }
    
    def design_allele_specific_primer(self, seq, pos, allele, probe_type, tail_length, target_tm):
        """设计等位基因特异性引物"""
        # 获取探针序列
        probe_seq = self.ksap_probes[probe_type]
        
        # 设计尾序列
        tail_seq = seq[max(0, pos-tail_length):pos]
        
        # 构建完整引物: 探针 + 尾序列 + 等位碱基
        primer_seq = probe_seq + tail_seq + allele
        
        # 计算Tm值
        primer_tm = mt.Tm_NN(Seq(primer_seq))
        
        return {
            'tail_primer': probe_seq + tail_seq,
            'allele_primer': tail_seq + allele,
            'sequence': primer_seq,
            'tm': round(primer_tm, 1)
        }
    
    def design_reverse_primer(self, seq, target_tm, tm_diff_threshold):
        """设计反向引物"""
        # 从序列末端开始寻找最佳引物
        best_primer = None
        best_diff = float('inf')
        
        min_len, max_len = self.config['primer_size_range']
        
        for length in range(min_len, max_len + 1):
            primer_seq = str(Seq(seq[-length:]).reverse_complement())
            primer_tm = mt.Tm_NN(Seq(primer_seq))
            tm_diff = abs(primer_tm - target_tm)
            
            if tm_diff < best_diff:
                best_diff = tm_diff
                best_primer = {
                    'sequence': primer_seq,
                    'tm': round(primer_tm, 1)
                }
                
                if best_diff <= tm_diff_threshold:  # 满足Tm偏差要求
                    break
        
        return best_primer if best_primer else {'sequence': '', 'tm': 0}
    
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
        """导出当前显示的结果"""
        if self.design_results is None or len(self.design_results) == 0:
            messagebox.showwarning("警告", "没有可导出的结果")
            return
            
        # 获取当前设计类型的名称用于默认文件名
        design_name = self.get_design_type_name()
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # 根据设计类型设置不同的默认文件名
        if self.current_design_type == "probe":
            default_name = f"探针设计结果_{timestamp}"
        elif self.current_design_type == "primer":
            default_name = f"扩增引物设计结果_{timestamp}"
        elif self.current_design_type == "kasp":
            default_name = f"KSAP引物设计结果_{timestamp}"
        else:
            default_name = f"分子设计结果_{timestamp}"
        
        # 弹出保存文件对话框
        path = filedialog.asksaveasfilename(
            initialdir=self.config['output_dir'],
            initialfile=default_name,
            defaultextension=".xlsx",
            filetypes=[
                ("Excel文件", "*.xlsx"),
                ("CSV文件", "*.csv"), 
                ("所有文件", "*.*")
            ],
            title="保存设计结果"
        )
        
        if not path:  # 用户取消了保存
            return
            
        try:
            # 根据文件扩展名选择保存格式
            if path.lower().endswith('.csv'):
                self.design_results.to_csv(path, index=False, encoding='utf_8_sig')  # 使用utf_8_sig编码避免中文乱码
            else:  # 默认为Excel格式
                self.design_results.to_excel(path, index=False)
                
            # 显示成功消息
            self.status_var.set(f"结果已成功保存到: {path}")
            messagebox.showinfo("导出成功", f"设计结果已保存到:\n{path}")
            
            # 尝试打开所在文件夹
            try:
                if os.name == 'nt':  # Windows系统
                    os.startfile(os.path.dirname(path))
                elif os.name == 'posix':  # Mac/Linux
                    os.system(f'open "{os.path.dirname(path)}"')
            except:
                pass
            
        except PermissionError:
            messagebox.showerror("导出失败", "文件可能正在被其他程序打开，请关闭后重试")
            self.status_var.set("导出失败 - 文件被占用")
        except Exception as e:
            messagebox.showerror("导出失败", f"保存文件时出错:\n{str(e)}")
            self.status_var.set("导出失败")


if __name__ == "__main__":
    root = tk.Tk()
    app = MolecularDesignTool(root)
    root.mainloop()

def create_frame(parent):
    frame = tk.Frame(parent)
    MolecularDesignTool(frame)
    return frame