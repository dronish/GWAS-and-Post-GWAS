# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, colorchooser
import pandas as pd
import numpy as np
import concurrent.futures
from scipy import stats

import gzip
import re

def parse_vcf(file_path):
    """自动适配简化和原始 VCF 格式的解析器"""
    try:
        opener = gzip.open if file_path.endswith('.gz') else open
        samples = []
        snp_data = {'CHROM': [], 'POS': [], 'ID': []}
        gt_index = 0

        with opener(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    samples = parts[9:]
                    for sample in samples:
                        snp_data[sample] = []
                    break

        if not samples:
            raise ValueError("VCF文件中未找到样本信息")

        with opener(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                var_id = parts[2] if parts[2] != '.' else f"{chrom}_{pos}"
                ref = parts[3]
                alts = parts[4].split(',')

                if len(ref) != 1 or any(len(alt) != 1 for alt in alts):
                    continue

                format_fields = parts[8].split(':') if len(parts) > 8 else []
                if 'GT' in format_fields:
                    gt_index = format_fields.index('GT')
                else:
                    gt_index = 0

                snp_data['CHROM'].append(chrom)
                snp_data['POS'].append(pos)
                snp_data['ID'].append(var_id)

                for i, sample_gt in enumerate(parts[9:]):
                    sample_name = samples[i]
                    gt_parts = sample_gt.split(':')
                    if gt_index >= len(gt_parts):
                        snp_data[sample_name].append(np.nan)
                        continue
                    gt = gt_parts[gt_index]
                    if gt in ('./.', '.', ''):
                        snp_data[sample_name].append(np.nan)
                    else:
                        alleles = re.split(r'[|/]', gt)
                        try:
                            allele1 = int(alleles[0])
                            allele2 = int(alleles[1]) if len(alleles) > 1 else allele1
                            dosage = allele1 + allele2
                            snp_data[sample_name].append(dosage)
                        except Exception:
                            snp_data[sample_name].append(np.nan)

        return pd.DataFrame(snp_data)

    except Exception as e:
        raise ValueError(f"解析VCF文件失败: {str(e)}")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import psutil
import os
from datetime import datetime
import webbrowser
import gc
import sys
from sklearn.decomposition import PCA
from scipy.linalg import lstsq
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

class GWASApp:
    def __init__(self, parent):
        self.root = parent 
        
        # 配置设置
        self.manhattan_colors = ['#1f77b4', '#ff7f0e']
        self.qq_line_color = 'red'
        self.sig_threshold = 5e-8
        self.suggestive_threshold = 1e-5
        self.point_size = 8
        self._threads = 4
        self._max_memory_gb = 8
        self.output_dir = os.path.expanduser("~")
        self._model_type = 'linear'
        self._pca_components = 3
        self.kinship_matrix = None
        
        # 数据存储
        self.genotype_data = None
        self.phenotype_data = None
        self.gwas_results = None
        self.current_figure = None
        self.vcf_path = None
        self.pca_result = None
        self.analysis_running = False

        # 初始化UI
        self.setup_ui()
        #self.setup_menu()

    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        try:
            self._threads = max(1, min(int(value), os.cpu_count() or 4))
            self.thread_var.set(str(self._threads))
        except ValueError:
            messagebox.showerror("错误", "请输入有效的线程数")

    @property
    def max_memory_gb(self):
        return self._max_memory_gb
    
    @max_memory_gb.setter
    def max_memory_gb(self, value):
        try:
            max_available = psutil.virtual_memory().total / (1024**3)
            self._max_memory_gb = max(1, min(float(value), max_available))
            self.mem_var.set(str(self._max_memory_gb))
        except ValueError:
            messagebox.showerror("错误", "请输入有效的内存限制")

    @property
    def model_type(self):
        return self._model_type
    
    @model_type.setter
    def model_type(self, value):
        if value in ['linear', 'mixed']:
            self._model_type = value
            self.model_var.set(value)

    @property
    def pca_components(self):
        return self._pca_components
    
    @pca_components.setter
    def pca_components(self, value):
        try:
            new_value = max(0, min(int(value), 10))
            if new_value != self._pca_components:
                self._pca_components = new_value
                self.pca_var.set(str(self._pca_components))
                # 当PCA组件数改变时，重置已有的PCA结果
                self.pca_result = None
        except ValueError:
            messagebox.showerror("错误", "请输入有效的PCA组件数")

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
        
        btn_load_pheno = tk.Button(control_frame, text="加载表型数据", command=self.load_phenotype)
        btn_load_pheno.pack(side=tk.LEFT, padx=5)
        
        btn_load_kinship = tk.Button(control_frame, text="加载亲缘矩阵", command=self.load_kinship_matrix)
        btn_load_kinship.pack(side=tk.LEFT, padx=5)
        
        # 分析控制按钮
        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(side=tk.LEFT, padx=10)
        
        self.btn_run = tk.Button(btn_frame, text="运行分析", command=self.run_gwas)
        self.btn_run.pack(side=tk.LEFT, padx=5)
        
        btn_stop = tk.Button(btn_frame, text="停止", command=self.stop_analysis)
        btn_stop.pack(side=tk.LEFT, padx=5)
        
        # 资源设置
        setting_frame = tk.Frame(control_frame)
        setting_frame.pack(side=tk.LEFT, padx=10)
        
        # 模型选择
        model_frame = tk.Frame(setting_frame)
        model_frame.pack(side=tk.LEFT, padx=5)
        tk.Label(model_frame, text="模型:").pack(side=tk.LEFT)
        self.model_var = tk.StringVar(value=self.model_type)
        self.model_menu = tk.OptionMenu(model_frame, self.model_var, 'linear', 'mixed', 
                                      command=lambda v: setattr(self, 'model_type', v))
        self.model_menu.pack(side=tk.LEFT)
        
        # PCA组件设置
        pca_frame = tk.Frame(setting_frame)
        pca_frame.pack(side=tk.LEFT, padx=5)
        tk.Label(pca_frame, text="PCA:").pack(side=tk.LEFT)
        self.pca_var = tk.StringVar(value=str(self.pca_components))
        pca_spin = tk.Spinbox(pca_frame, from_=0, to=10, textvariable=self.pca_var, width=3,
                             command=lambda: setattr(self, 'pca_components', self.pca_var.get()))
        pca_spin.pack(side=tk.LEFT)
        
        # 线程设置
        thread_frame = tk.Frame(setting_frame)
        thread_frame.pack(side=tk.LEFT, padx=5)
        tk.Label(thread_frame, text="线程数:").pack(side=tk.LEFT)
        self.thread_var = tk.StringVar(value=str(self.threads))
        thread_entry = tk.Entry(thread_frame, textvariable=self.thread_var, width=3)
        thread_entry.pack(side=tk.LEFT)
        thread_entry.bind("<FocusOut>", lambda e: setattr(self, 'threads', self.thread_var.get()))
        
        # 内存设置
        mem_frame = tk.Frame(setting_frame)
        mem_frame.pack(side=tk.LEFT, padx=5)
        tk.Label(mem_frame, text="最大内存(GB):").pack(side=tk.LEFT)
        self.mem_var = tk.StringVar(value=str(self.max_memory_gb))
        mem_entry = tk.Entry(mem_frame, textvariable=self.mem_var, width=3)
        mem_entry.pack(side=tk.LEFT)
        mem_entry.bind("<FocusOut>", lambda e: setattr(self, 'max_memory_gb', self.mem_var.get()))
        
        # 进度条
        progress_frame = tk.Frame(main_frame)
        progress_frame.pack(fill=tk.X, pady=5)
        
        self.progress = ttk.Progressbar(progress_frame, orient=tk.HORIZONTAL, length=400, mode='determinate')
        self.progress.pack(fill=tk.X)
        self.progress_label = tk.Label(progress_frame, text="准备就绪")
        self.progress_label.pack()
        
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
        
        # 图形+按钮区域整体容器
        viz_frame_container = tk.Frame(main_frame)
        viz_frame_container.pack(fill=tk.BOTH, expand=True)

        # 右上角按钮区域
        viz_button_frame = tk.Frame(viz_frame_container)
        viz_button_frame.pack(anchor='ne', pady=2)
        tk.Button(viz_button_frame, text="曼哈顿图", command=self.plot_manhattan).pack(side=tk.LEFT, padx=2)
        tk.Button(viz_button_frame, text="QQ图", command=self.plot_qq).pack(side=tk.LEFT, padx=2)
        tk.Button(viz_button_frame, text="绘图设置", command=self.open_plot_settings).pack(side=tk.LEFT, padx=2)
        tk.Button(viz_button_frame, text="保存结果", command=self.save_results).pack(side=tk.LEFT, padx=2)

        # 图形显示区域
        self.viz_frame = tk.Frame(viz_frame_container, bg='white')
        self.viz_frame.pack(fill=tk.BOTH, expand=True)

    def setup_menu(self):
        menubar = tk.Menu(self.root)
        
        # 文件菜单
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="加载VCF", command=self.load_vcf)
        file_menu.add_command(label="加载表型", command=self.load_phenotype)
        file_menu.add_command(label="加载亲缘矩阵", command=self.load_kinship_matrix)
        file_menu.add_separator()
        file_menu.add_command(label="保存结果", command=self.save_results)
        file_menu.add_separator()
        file_menu.add_command(label="退出", command=self.root.quit)
        menubar.add_cascade(label="文件", menu=file_menu)
        
        # 分析菜单
        analysis_menu = tk.Menu(menubar, tearoff=0)
        analysis_menu.add_command(label="运行分析", command=self.run_gwas)
        analysis_menu.add_command(label="停止分析", command=self.stop_analysis)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="计算PCA", command=self.calculate_pca)
        analysis_menu.add_command(label="计算亲缘矩阵", command=self.calculate_kinship)
        menubar.add_cascade(label="分析", menu=analysis_menu)
        
        # 设置菜单
        settings_menu = tk.Menu(menubar, tearoff=0)
        settings_menu.add_command(label="绘图设置", command=self.open_plot_settings)
        settings_menu.add_command(label="模型设置", command=self.open_model_settings)
        settings_menu.add_command(label="输出目录", command=self.set_output_dir)
        menubar.add_cascade(label="设置", menu=settings_menu)
        
        # 帮助菜单
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="使用帮助", command=self.show_help)
        help_menu.add_command(label="关于", command=self.show_about)
        menubar.add_cascade(label="帮助", menu=help_menu)
        
        self.root.config(menu=menubar)

    
    def load_vcf(self):
        if self.analysis_running:
            messagebox.showwarning("警告", "请先停止当前分析")
            return

        if not self.check_memory_usage():
            return

        path = filedialog.askopenfilename(
            initialdir=self.output_dir,
            filetypes=[("VCF文件", "*.vcf *.vcf.gz"), ("所有文件", "*.*")])

        if not path:
            return

        self.vcf_path = path
        try:
            self.progress_label.config(text="正在加载VCF文件...")
            self.progress["value"] = 0
            self.root.update()

            self.genotype_data = parse_vcf(path)
            self.kinship_matrix = None
            self.pca_result = None
            gc.collect()

            missing_rate = self.genotype_data.iloc[:, 3:].isna().mean(axis=1)
            self.genotype_data = self.genotype_data[missing_rate <= 0.2]

            mem_usage = self.genotype_data.memory_usage(deep=True).sum() / (1024 ** 3)
            self.progress["value"] = 100
            self.progress_label.config(text=f"VCF加载完成 (占用内存: {mem_usage:.2f}GB)")
            messagebox.showinfo("成功", f"成功加载 {len(self.genotype_data)} 个SNP\n占用内存: {mem_usage:.2f}GB")

        except Exception as e:
            messagebox.showerror("错误", f"加载VCF失败:\n{str(e)}")
            self.progress_label.config(text="加载失败")


    def load_phenotype(self):
        if self.analysis_running:
            messagebox.showwarning("警告", "请先停止当前分析")
            return
            
        path = filedialog.askopenfilename(
            initialdir=self.output_dir,
            filetypes=[("CSV文件", "*.csv"), ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            self.progress_label.config(text="正在加载表型数据...")
            self.root.update()
            
            if path.endswith('.csv'):
                df = pd.read_csv(path, index_col=0)
            else:
                df = pd.read_excel(path, index_col=0)
            
            df = df.dropna(how='all').dropna(axis=1, how='all')
            self.phenotype_data = df
            
            self.progress_label.config(text="表型数据加载完成")
            messagebox.showinfo("成功", f"加载 {len(self.phenotype_data)} 个样本的表型数据")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载表型数据失败:\n{str(e)}")
            self.progress_label.config(text="加载失败")

    def load_kinship_matrix(self):
        path = filedialog.askopenfilename(
            initialdir=self.output_dir,
            filetypes=[("CSV文件", "*.csv"), ("文本文件", "*.txt"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            if path.endswith('.csv'):
                kinship = pd.read_csv(path, index_col=0)
            else:
                kinship = pd.read_csv(path, sep='\t', index_col=0)
            
            if not np.allclose(kinship.values, kinship.values.T):
                messagebox.showwarning("警告", "亲缘矩阵不对称，将自动对称化")
                kinship = (kinship + kinship.T) / 2
            
            # 检查矩阵是否有效
            if np.any(np.isnan(kinship.values)) or np.any(np.isinf(kinship.values)):
                messagebox.showerror("错误", "亲缘矩阵包含无效值(NaN或inf)")
                return
                
            self.kinship_matrix = kinship
            messagebox.showinfo("成功", f"成功加载亲缘矩阵\n样本数: {len(kinship)}")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载亲缘矩阵失败:\n{str(e)}")

    def calculate_pca(self):
        if self.genotype_data is None:
            messagebox.showwarning("警告", "请先加载VCF文件")
            return
            
        try:
            if self.pca_components <= 0:
                messagebox.showinfo("信息", "PCA组件数设置为0，将不计算PCA")
                self.pca_result = None
                return
                
            samples = [col for col in self.genotype_data.columns if col not in ['CHROM', 'POS', 'ID']]
            geno_matrix = self.genotype_data[samples].values.T
            
            # 填充缺失值 - 按列均值填充
            col_means = np.nanmean(geno_matrix, axis=0)
            geno_matrix = np.where(np.isnan(geno_matrix), col_means, geno_matrix)
            
            # 检查是否还有缺失值
            if np.any(np.isnan(geno_matrix)):
                messagebox.showwarning("警告", "基因型数据中存在过多缺失值")
                return
                
            # 标准化
            geno_std = (geno_matrix - np.mean(geno_matrix, axis=0)) / np.std(geno_matrix, axis=0)
            geno_std = np.nan_to_num(geno_std, nan=0.0, posinf=0.0, neginf=0.0)
            
            # 计算PCA
            pca = PCA(n_components=self.pca_components)
            self.pca_result = pca.fit_transform(geno_std)
            
            messagebox.showinfo("成功", f"PCA计算完成\n解释方差: {np.sum(pca.explained_variance_ratio_):.2f}")
            
        except Exception as e:
            messagebox.showerror("错误", f"计算PCA失败:\n{str(e)}")

    def calculate_kinship(self):
        if self.genotype_data is None:
            messagebox.showwarning("警告", "请先加载VCF文件")
            return False
            
        try:
            self.progress_label.config(text="正在计算亲缘矩阵...")
            self.progress["value"] = 0
            self.root.update()
            
            samples = [col for col in self.genotype_data.columns if col not in ['CHROM', 'POS', 'ID']]
            geno_matrix = self.genotype_data[samples].values.T
            
            # 填充缺失值 - 按列均值填充
            col_means = np.nanmean(geno_matrix, axis=0)
            geno_matrix = np.where(np.isnan(geno_matrix), col_means, geno_matrix)
            
            # 检查是否还有缺失值
            if np.any(np.isnan(geno_matrix)):
                messagebox.showwarning("警告", "基因型数据中存在过多缺失值")
                return False
                
            # 标准化
            geno_std = (geno_matrix - np.mean(geno_matrix, axis=0)) / np.std(geno_matrix, axis=0)
            geno_std = np.nan_to_num(geno_std, nan=0.0, posinf=0.0, neginf=0.0)
            
            self.progress["value"] = 50
            self.progress_label.config(text="正在计算基因组关系矩阵...")
            self.root.update()
            
            # 计算GRM (VanRaden方法)
            grm = np.dot(geno_std, geno_std.T) / geno_std.shape[1]
            
            # 对角元素调整为1
            np.fill_diagonal(grm, 1)
            
            # 检查GRM是否有效
            if np.any(np.isnan(grm)) or np.any(np.isinf(grm)):
                messagebox.showwarning("警告", "亲缘矩阵计算产生无效值")
                return False
                
            self.kinship_matrix = pd.DataFrame(grm, index=samples, columns=samples)
            
            self.progress["value"] = 100
            self.progress_label.config(text="亲缘矩阵计算完成")
            messagebox.showinfo("成功", "亲缘矩阵计算完成")
            return True
            
        except Exception as e:
            messagebox.showerror("错误", f"计算亲缘矩阵失败:\n{str(e)}")
            self.progress_label.config(text="计算失败")
            return False

    def handle_missing_kinship(self):
        dialog = tk.Toplevel(self.root)
        dialog.title("缺少亲缘矩阵")
        dialog.geometry("400x200")
        
        tk.Label(dialog, text="您选择了混合线性模型但未提供亲缘矩阵", font=('Arial', 12)).pack(pady=10)
        
        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=20)
        
        result = {"choice": None}
        
        def set_choice(choice):
            result["choice"] = choice
            dialog.destroy()
        
        tk.Button(btn_frame, text="使用简单线性模型", 
                 command=lambda: set_choice("linear")).pack(side=tk.LEFT, padx=10)
        
        tk.Button(btn_frame, text="自动计算亲缘矩阵", 
                 command=lambda: set_choice("calculate")).pack(side=tk.LEFT, padx=10)
        
        tk.Button(btn_frame, text="取消分析", 
                 command=lambda: set_choice("cancel")).pack(side=tk.LEFT, padx=10)
        
        dialog.transient(self.root)
        dialog.grab_set()
        self.root.wait_window(dialog)
        
        return result["choice"]

    def run_gwas(self):
        if self.analysis_running:
            return
            
        # 确保使用当前设置
        self.threads = self.thread_var.get()
        self.max_memory_gb = self.mem_var.get()
        self.model_type = self.model_var.get()
        self.pca_components = self.pca_var.get()
        
        if not self.check_memory_usage():
            return
            
        if self.genotype_data is None:
            messagebox.showwarning("警告", "请先加载VCF文件")
            return
            
        if self.phenotype_data is None:
            messagebox.showwarning("警告", "请先加载表型数据")
            return
            
        # 处理混合模型但没有亲缘矩阵的情况
        if self.model_type == 'mixed' and self.kinship_matrix is None:
            choice = self.handle_missing_kinship()
            
            if choice == "cancel":
                return
            elif choice == "linear":
                self.model_type = 'linear'
            elif choice == "calculate":
                if not self.calculate_kinship():
                    return
        
        try:
            self.analysis_running = True
            self.btn_run.config(state=tk.DISABLED)
            self.progress_label.config(text="正在准备分析数据...")
            self.progress["value"] = 0
            self.root.update()
            
            # 获取样本列表
            samples = [col for col in self.genotype_data.columns if col not in ['CHROM', 'POS', 'ID']]
            
            # 确保表型数据是Series且有名称为"pheno"
            if isinstance(self.phenotype_data, pd.DataFrame):
                y = self.phenotype_data.iloc[:, 0].rename("pheno")
            else:
                y = self.phenotype_data.copy()
                y.name = "pheno"
            
            # 找出共同样本
            common_samples = list(set(samples) & set(y.index))
            
            if len(common_samples) < 10:
                messagebox.showerror("错误", f"共同样本数不足(仅{len(common_samples)}个)")
                self.analysis_running = False
                self.btn_run.config(state=tk.NORMAL)
                return
                
            # 如果PCA组件数>0且没有PCA结果或PCA组件数不匹配，则重新计算PCA
            if self.pca_components > 0 and (self.pca_result is None or self.pca_result.shape[1] != self.pca_components):
                self.calculate_pca()
                
            # 准备协变量（PCA）
            covariates = None
            if self.pca_components > 0 and self.pca_result is not None:
                # 确保PCA样本顺序与共同样本一致
                sample_indices = [samples.index(s) for s in common_samples]
                covariates = self.pca_result[sample_indices, :self.pca_components]
                # 标准化PCA成分
                covariates = (covariates - np.mean(covariates, axis=0)) / np.std(covariates, axis=0)
                covariates = np.nan_to_num(covariates, nan=0.0, posinf=0.0, neginf=0.0)
            
            # 准备亲缘矩阵（混合模型）
            kinship = None
            if self.model_type == 'mixed' and self.kinship_matrix is not None:
                # 确保亲缘矩阵样本与共同样本匹配
                kinship_samples = [s for s in common_samples if s in self.kinship_matrix.index]
                if len(kinship_samples) < 10:
                    messagebox.showerror("错误", f"亲缘矩阵样本匹配不足(仅{len(kinship_samples)}个)")
                    self.analysis_running = False
                    self.btn_run.config(state=tk.NORMAL)
                    return
                
                # 确保表型数据与亲缘矩阵样本一致
                y = y[kinship_samples]
                kinship = self.kinship_matrix.loc[kinship_samples, kinship_samples].values
                
                # 更新共同样本列表
                common_samples = kinship_samples
            
            # 多线程分析
            results = []
            total_snps = len(self.genotype_data)
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
                futures = []
                for i, row in self.genotype_data.iterrows():
                    if not self.analysis_running:
                        break
                        
                    futures.append(executor.submit(
                        self.process_snp, 
                        row, 
                        common_samples, 
                        y,
                        covariates,
                        kinship
                    ))
                    
                    if i % 500 == 0:
                        self.progress["value"] = (i / total_snps) * 90
                        self.progress_label.config(text=f"分析进度: {i}/{total_snps}")
                        self.root.update()
                
                for j, future in enumerate(concurrent.futures.as_completed(futures)):
                    if not self.analysis_running:
                        break
                        
                    res = future.result()
                    if res:
                        results.append(res)
                        
                    if j % 500 == 0:
                        self.progress["value"] = 90 + (j / len(futures)) * 10
                        self.progress_label.config(text=f"整理结果: {j}/{len(futures)}")
                        self.root.update()
            
            if results and self.analysis_running:
                self.gwas_results = pd.DataFrame(results).sort_values(by='P')
                self.plot_manhattan()
                
                # 计算lambda GC
                observed_pvals = self.gwas_results['P'].dropna()
                if len(observed_pvals) > 0:
                    chi2 = stats.chi2.ppf(1 - observed_pvals, 1)
                    lambda_gc = np.median(chi2) / 0.4549
                    messagebox.showinfo("分析完成", 
                                      f"分析完成，共发现 {len(results)} 个有效SNP\n"
                                      f"Lambda GC: {lambda_gc:.4f}")
                else:
                    messagebox.showinfo("分析完成", "分析完成，但未发现显著SNP")
                
                self.progress["value"] = 100
                self.progress_label.config(text="分析完成")
            else:
                self.progress_label.config(text="分析已停止")
                
            self.update_result_table()
            
        except Exception as e:
            messagebox.showerror("错误", f"分析失败:\n{str(e)}")
            self.progress_label.config(text="分析失败")
            
        finally:
            self.analysis_running = False
            self.btn_run.config(state=tk.NORMAL)

    def process_snp(self, row, samples, y, covariates=None, kinship=None):
        try:
            # 获取基因型数据 - 确保转换为数值型numpy数组
            geno = pd.to_numeric(row[samples], errors='coerce').to_numpy()
            
            # 检查基因型数据是否有效
            if np.all(np.isnan(geno)):
                return None
                
            # 创建数据框并删除缺失值
            data = pd.DataFrame({'geno': geno, 'pheno': y}).dropna()
            
            if len(data) < 10 or data['geno'].nunique() < 2:
                return None
                
            # 标准化基因型 - 确保转换为float类型
            geno_mean = data['geno'].mean()
            geno_std = data['geno'].std()
            if geno_std == 0:  # 防止除以0
                return None
                
            geno_std = ((data['geno'] - geno_mean) / geno_std).astype(np.float64)
            geno_std = np.nan_to_num(geno_std, nan=0.0, posinf=0.0, neginf=0.0)
            
            # 检查标准化后的基因型是否有效
            if np.any(np.isnan(geno_std)) or np.any(np.isinf(geno_std)):
                return None
            
            # 准备设计矩阵
            X = geno_std.reshape(-1, 1)  # 确保是float类型的numpy数组
            
            # 添加协变量（PCA）
            if covariates is not None:
                # 获取当前样本在共同样本中的索引
                sample_indices = [samples.index(s) for s in data.index]
                if len(sample_indices) != len(data):
                    return None
                covar = covariates[sample_indices].astype(np.float64)
                
                # 检查协变量是否有效
                if np.any(np.isnan(covar)) or np.any(np.isinf(covar)):
                    return None
                    
                # 添加截距项和PCA成分
                X = np.hstack([X, np.ones((len(data), 1), dtype=np.float64), covar])
            else:
                # 仅添加截距项
                X = np.hstack([X, np.ones((len(data), 1), dtype=np.float64)])
            
            # 检查设计矩阵是否有效
            if np.any(np.isnan(X)) or np.any(np.isinf(X)):
                return None
            
            # 中心化表型 - 确保使用float类型
            y_centered = (data['pheno'].to_numpy() - data['pheno'].mean()).astype(np.float64)
            y_centered = np.nan_to_num(y_centered, nan=0.0, posinf=0.0, neginf=0.0)
            
            # 检查表型是否有效
            if np.any(np.isnan(y_centered)) or np.any(np.isinf(y_centered)):
                return None
            
            # 选择模型
            if kinship is not None:
                # 混合线性模型
                sample_indices = [samples.index(s) for s in data.index]
                current_kinship = kinship[np.ix_(sample_indices, sample_indices)].astype(np.float64)
                current_kinship = np.nan_to_num(current_kinship, nan=0.0, posinf=0.0, neginf=0.0)
                
                # 检查亲缘矩阵是否有效
                if np.any(np.isnan(current_kinship)) or np.any(np.isinf(current_kinship)):
                    return None
                    
                result = self.run_mixed_model(X, y_centered, current_kinship)
            else:
                # 简单线性模型
                result = self.run_linear_model(X, y_centered)
        
            if result is None:
                return None
                
            return {
                'CHR': row['CHROM'],
                'POS': row['POS'],
                'SNP': row['ID'],
                'Beta': float(result['beta']),
                'SE': float(result['se']),
                'P': float(result['p_value'])
            }
            
        except Exception as e:
            print(f"处理SNP {row['ID']} 出错: {str(e)}")
            return None

    def run_linear_model(self, X, y):
        try:
            model = OLS(y, X)
            results = model.fit()
            
            return {
                'beta': results.params[0],  # 第一个系数是SNP效应
                'se': results.bse[0],
                'p_value': results.pvalues[0]
            }
        except Exception as e:
            print(f"线性模型计算错误: {str(e)}")
            return None

    def run_mixed_model(self, X, y, K):
        """改进的混合线性模型实现"""
        try:
            n = len(y)
            
            # 1. 对K进行谱分解
            S, U = np.linalg.eigh(K)
            S = np.maximum(S, 1e-6)  # 确保正定性
            
            # 2. 转换变量
            y_t = U.T @ y
            X_t = U.T @ X
            
            # 3. 使用EM算法估计方差组分
            # 初始值
            sigma_g = np.var(y) * 0.5
            sigma_e = np.var(y) * 0.5
            prev_loglik = -np.inf
            
            for _ in range(20):  # 最多20次迭代
                # 构建权重矩阵
                weights = 1.0 / (S * sigma_g + sigma_e)
                W = np.diag(weights)
                
                # 计算beta估计
                XWX = X_t.T @ W @ X_t
                XWy = X_t.T @ W @ y_t
                beta = lstsq(XWX, XWy)[0]
                
                # 计算残差和更新方差组分
                resid = y_t - X_t @ beta
                sigma_g = max(1e-6, (resid.T @ np.diag(S * weights**2) @ resid) / np.sum(S * weights))
                sigma_e = max(1e-6, (resid.T @ np.diag(weights**2) @ resid) / np.sum(weights))
                
                # 检查收敛
                loglik = -0.5 * (np.sum(np.log(S * sigma_g + sigma_e)) + np.sum(weights * resid**2))
                if abs(loglik - prev_loglik) < 1e-4:
                    break
                prev_loglik = loglik
            
            # 最终beta估计和标准误
            weights = 1.0 / (S * sigma_g + sigma_e)
            W = np.diag(weights)
            XWX = X_t.T @ W @ X_t
            XWy = X_t.T @ W @ y_t
            beta = lstsq(XWX, XWy)[0]
            
            # 计算协方差矩阵
            cov_beta = np.linalg.inv(XWX)
            se = np.sqrt(np.diag(cov_beta))[0]  # 第一个是SNP效应
            
            # 计算p值
            t_stat = beta[0] / se
            p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-X_t.shape[1]))
            
            return {
                'beta': beta[0],
                'se': se,
                'p_value': p_value,
                'sigma_g': sigma_g,
                'sigma_e': sigma_e
            }
            
        except Exception as e:
            print(f"混合模型计算错误: {str(e)}")
            return None

    def plot_manhattan(self):
        if self.gwas_results is None or self.gwas_results.empty:
            messagebox.showwarning("警告", "没有结果可绘制")
            return
            
        try:
            df = self.gwas_results.copy()
            df['-log10(P)'] = -np.log10(df['P'])
            df['CHR'] = df['CHR'].astype(str)
            
            def chr_key(x):
                try:
                    return int(x.replace('chr', '').replace('Chr', ''))
                except ValueError:
                    return float('inf')
            
            df = df.sort_values('CHR', key=lambda x: x.map(chr_key))
            
            chr_groups = df.groupby('CHR', observed=True, sort=False)
            chr_lengths = chr_groups['POS'].max()
            chr_offset = chr_lengths.cumsum().shift(1).fillna(0)
            df['abs_pos'] = df['POS'] + df['CHR'].map(chr_offset)
            
            fig, ax = plt.subplots(figsize=(12, 6))
            colors = self.manhattan_colors
            
            for i, (chr_name, group) in enumerate(chr_groups):
                color = colors[i % len(colors)]
                ax.scatter(group['abs_pos'], group['-log10(P)'], 
                          color=color, s=self.point_size)
            
            sig_level = -np.log10(self.sig_threshold)
            sug_level = -np.log10(self.suggestive_threshold)
            ax.axhline(sig_level, color='r', linestyle='--', linewidth=1)
            ax.axhline(sug_level, color='g', linestyle='--', linewidth=1)
            
            chr_pos = df.groupby('CHR', observed=True, sort=False)['abs_pos'].median()
            sorted_chr = sorted(chr_pos.index, key=chr_key)
            ax.set_xticks([chr_pos[chr] for chr in sorted_chr])
            ax.set_xticklabels(sorted_chr)
            
            ax.set_xlabel("chr")
            ax.set_ylabel("-log10(pvalue)")
            ax.set_title(f"Manhattan plot ({self.model_type}model)")
            
            plt.tight_layout()
            self.show_plot(fig)
            
        except Exception as e:
            messagebox.showerror("错误", f"绘制曼哈顿图失败:\n{str(e)}")

    def plot_qq(self):
        if self.gwas_results is None or self.gwas_results.empty:
            messagebox.showwarning("警告", "没有结果可绘制")
            return
            
        try:
            pvals = self.gwas_results['P'].dropna().sort_values()
            expected = -np.log10(np.linspace(1/len(pvals), 1, len(pvals)))
            observed = -np.log10(pvals)
            
            # 计算lambda GC
            chi2 = stats.chi2.ppf(1 - pvals, 1)
            lambda_gc = np.median(chi2) / 0.4549
            
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(expected, observed, s=self.point_size)
            
            ax.plot([0, max(expected)], [0, max(expected)], 
                    color=self.qq_line_color,
                    linestyle='--')
            
            ax.set_xlabel("Expected -log10(P)")
            ax.set_ylabel("Observed -log10(P)")
            ax.set_title(f"QQplot ({self.model_type}model)\nλGC={lambda_gc:.4f}")
            
            plt.tight_layout()
            self.show_plot(fig)
            
        except Exception as e:
            messagebox.showerror("错误", f"绘制QQ图失败:\n{str(e)}")

    def show_plot(self, fig):
        for widget in self.viz_frame.winfo_children():
            widget.destroy()
        
        canvas = FigureCanvasTkAgg(fig, master=self.viz_frame)
        canvas.draw()
        
        toolbar = NavigationToolbar2Tk(canvas, self.viz_frame)
        toolbar.update()
        
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.current_figure = fig

    def open_plot_settings(self):
        settings_win = tk.Toplevel(self.root)
        settings_win.title("绘图设置")
        settings_win.geometry("400x300")
        
        tk.Label(settings_win, text="显著性阈值:").grid(row=0, column=0, padx=10, pady=5, sticky='e')
        sig_entry = tk.Entry(settings_win)
        sig_entry.insert(0, str(self.sig_threshold))
        sig_entry.grid(row=0, column=1, padx=10, pady=5, sticky='w')
        
        tk.Label(settings_win, text="提示性阈值:").grid(row=1, column=0, padx=10, pady=5, sticky='e')
        sug_entry = tk.Entry(settings_win)
        sug_entry.insert(0, str(self.suggestive_threshold))
        sug_entry.grid(row=1, column=1, padx=10, pady=5, sticky='w')
        
        tk.Label(settings_win, text="点大小:").grid(row=2, column=0, padx=10, pady=5, sticky='e')
        size_entry = tk.Entry(settings_win)
        size_entry.insert(0, str(self.point_size))
        size_entry.grid(row=2, column=1, padx=10, pady=5, sticky='w')
        
        tk.Label(settings_win, text="曼哈顿图颜色1:").grid(row=3, column=0, padx=10, pady=5, sticky='e')
        color1_btn = tk.Button(settings_win, bg=self.manhattan_colors[0], width=10,
                              command=lambda: self.choose_color(color1_btn, 'manhattan_colors', 0))
        color1_btn.grid(row=3, column=1, padx=10, pady=5, sticky='w')
        
        tk.Label(settings_win, text="曼哈顿图颜色2:").grid(row=4, column=0, padx=10, pady=5, sticky='e')
        color2_btn = tk.Button(settings_win, bg=self.manhattan_colors[1], width=10,
                              command=lambda: self.choose_color(color2_btn, 'manhattan_colors', 1))
        color2_btn.grid(row=4, column=1, padx=10, pady=5, sticky='w')
        
        tk.Label(settings_win, text="QQ图参考线颜色:").grid(row=5, column=0, padx=10, pady=5, sticky='e')
        qq_color_btn = tk.Button(settings_win, bg=self.qq_line_color, width=10,
                                command=lambda: self.choose_color(qq_color_btn, 'qq_line_color'))
        qq_color_btn.grid(row=5, column=1, padx=10, pady=5, sticky='w')
        
        save_btn = tk.Button(settings_win, text="保存设置", 
                            command=lambda: self.save_plot_settings(
                                sig_entry.get(),
                                sug_entry.get(),
                                size_entry.get(),
                                settings_win))
        save_btn.grid(row=6, column=0, columnspan=2, pady=10)

    def save_plot_settings(self, sig_thresh, sug_thresh, point_size, win):
        try:
            self.sig_threshold = float(sig_thresh)
            self.suggestive_threshold = float(sug_thresh)
            self.point_size = int(point_size)
            messagebox.showinfo("成功", "绘图设置已保存")
            win.destroy()
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数值")

    def open_model_settings(self):
        settings_win = tk.Toplevel(self.root)
        settings_win.title("模型设置")
        settings_win.geometry("300x200")
        
        tk.Label(settings_win, text="分析模型:").grid(row=0, column=0, padx=10, pady=10, sticky='e')
        model_var = tk.StringVar(value=self.model_type)
        model_menu = tk.OptionMenu(settings_win, model_var, 'linear', 'mixed',
                                 command=lambda v: setattr(self, 'model_type', v))
        model_menu.grid(row=0, column=1, padx=10, pady=10, sticky='w')
        
        tk.Label(settings_win, text="PCA组件数:").grid(row=1, column=0, padx=10, pady=10, sticky='e')
        pca_var = tk.StringVar(value=str(self.pca_components))
        pca_spin = tk.Spinbox(settings_win, from_=0, to=10, textvariable=pca_var,
                             command=lambda: setattr(self, 'pca_components', pca_var.get()))
        pca_spin.grid(row=1, column=1, padx=10, pady=10, sticky='w')
        
        save_btn = tk.Button(settings_win, text="保存设置", 
                            command=lambda: settings_win.destroy())
        save_btn.grid(row=2, column=0, columnspan=2, pady=10)

    def choose_color(self, button, config_key, index=None):
        color = colorchooser.askcolor(title="选择颜色")[1]
        if color:
            button.config(bg=color)
            if index is not None:
                if config_key == 'manhattan_colors':
                    self.manhattan_colors[index] = color
            else:
                if config_key == 'qq_line_color':
                    self.qq_line_color = color

    def update_result_table(self):
        if self.gwas_results is None or self.gwas_results.empty:
            return
            
        self.result_tree.delete(*self.result_tree.get_children())
        self.result_tree["columns"] = list(self.gwas_results.columns)
        self.result_tree["show"] = "headings"
        
        col_widths = {'CHR': 60, 'POS': 80, 'SNP': 120, 'Beta': 80, 'SE': 80, 'P': 100}
        for col in self.gwas_results.columns:
            self.result_tree.heading(col, text=col)
            self.result_tree.column(col, width=col_widths.get(col, 100), anchor=tk.CENTER)
        
        for _, row in self.gwas_results.iterrows():
            self.result_tree.insert("", tk.END, values=list(row))

    def save_results(self):
        if self.gwas_results is None or self.gwas_results.empty:
            messagebox.showwarning("警告", "没有结果可保存")
            return
            
        default_name = f"gwas_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        path = filedialog.asksaveasfilename(
            initialdir=self.output_dir,
            initialfile=default_name,
            defaultextension=".csv",
            filetypes=[("CSV文件", "*.csv"), ("Excel文件", "*.xlsx"), ("所有文件", "*.*")])
        
        if not path:
            return
            
        try:
            if path.endswith('.csv'):
                self.gwas_results.to_csv(path, index=False)
            else:
                self.gwas_results.to_excel(path, index=False)
                
            messagebox.showinfo("成功", f"结果已保存到:\n{path}")
        except Exception as e:
            messagebox.showerror("错误", f"保存失败:\n{str(e)}")

    def set_output_dir(self):
        dir_path = filedialog.askdirectory(initialdir=self.output_dir)
        if dir_path:
            self.output_dir = dir_path
            messagebox.showinfo("成功", f"输出目录已设置为:\n{dir_path}")

    def stop_analysis(self):
        if self.analysis_running:
            self.analysis_running = False
            self.progress_label.config(text="分析已停止")
            messagebox.showinfo("信息", "正在停止分析...")
        else:
            messagebox.showinfo("信息", "没有正在运行的分析")

    def check_memory_usage(self):
        mem = psutil.virtual_memory()
        available_gb = mem.available / (1024 ** 3)
        
        if available_gb < self.max_memory_gb * 0.8:
            msg = (f"可用内存不足！\n"
                  f"当前可用: {available_gb:.1f}GB\n"
                  f"设置限制: {self.max_memory_gb}GB\n"
                  "请关闭其他程序或增加内存限制")
            messagebox.showwarning("内存警告", msg)
            return False
        return True

    def show_about(self):
        about_text = """GWAS分析工具 v3.5

功能特点:
- 支持VCF和表型数据加载
- 支持简单线性模型和混合线性模型
- 支持PCA作为协变量
- 自动计算亲缘关系矩阵
- 多线程GWAS分析
- 曼哈顿图和QQ图可视化

作者: dronish
联系方式: q15215390195@163.com"""
        messagebox.showinfo("关于", about_text)

    def show_help(self):
        help_url = "https://github.com/dronish/GWAS-and-Post-GWAS"
        webbrowser.open(help_url)

if __name__ == "__main__":
    root = tk.Tk()
    app = GWASApp(root)
    root.mainloop()

def create_frame(parent):
    frame = tk.Frame(parent)
    GWASApp(frame)  # 实例化 GUI 应用到这个 Frame
    return frame
