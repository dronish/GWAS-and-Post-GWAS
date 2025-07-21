# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk, messagebox
import webbrowser
import os
import subprocess
import sys
from PIL import Image, ImageTk
import gwas
import primer
import snpann

SOFTWARE_VERSION = "v1.0"
AUTHOR = "dronish"
CONTACT = "q15215390195@163.com"
Github = "https://github.com/dronish/GWAS-and-Post-GWAS"

def show_about():
    messagebox.showinfo("软件信息",
        f"优异等位基因快速提取工具\n版本：{SOFTWARE_VERSION}\n作者：{AUTHOR}\n联系方式：{CONTACT}\nGithub: {Github}"
    )

def open_folder(path):
    if sys.platform.startswith('darwin'):
        subprocess.call(['open', path])
    elif os.name == 'nt':
        os.startfile(path)
    elif os.name == 'posix':
        subprocess.call(['xdg-open', path])
    else:
        messagebox.showerror("错误", "当前系统不支持自动打开文件夹")

def open_instructions():
    path = os.path.abspath("说明")
    if os.path.exists(path):
        open_folder(path)
    else:
        messagebox.showerror("错误", "说明文件夹不存在")

def open_examples():
    path = os.path.abspath("示例")
    if os.path.exists(path):
        open_folder(path)
    else:
        messagebox.showerror("错误", "示例文件夹不存在")

def open_primer_blast():
    webbrowser.open("https://www.ncbi.nlm.nih.gov/tools/primer-blast/")

def main():
    root = tk.Tk()
    root.title("优异等位基因快速提取工具")
    root.geometry("1200x900")

    menubar = tk.Menu(root)

    welcome_menu = tk.Menu(menubar, tearoff=0)
    welcome_menu.add_command(label="欢迎页", command=lambda: notebook.select(welcome_tab))
    menubar.add_cascade(label="欢迎", menu=welcome_menu)

    analysis_menu = tk.Menu(menubar, tearoff=0)
    analysis_menu.add_command(label="GWAS分析", command=lambda: notebook.select(gwas_tab))
    analysis_menu.add_command(label="SNP注释", command=lambda: notebook.select(snp_tab))
    analysis_menu.add_command(label="引物/探针设计", command=lambda: notebook.select(primer_tab))
    menubar.add_cascade(label="分析", menu=analysis_menu)

    helpmenu = tk.Menu(menubar, tearoff=0)
    helpmenu.add_command(label="软件信息", command=show_about)
    helpmenu.add_command(label="使用说明", command=open_instructions)
    helpmenu.add_command(label="示例文件", command=open_examples)
    menubar.add_cascade(label="帮助", menu=helpmenu)

    blastmenu = tk.Menu(menubar, tearoff=0)
    blastmenu.add_command(label="打开NCBI引物BLAST", command=open_primer_blast)
    menubar.add_cascade(label="BLAST", menu=blastmenu)

    root.config(menu=menubar)

    global notebook, welcome_tab, gwas_tab, snp_tab, primer_tab
    notebook = ttk.Notebook(root)
    notebook.pack(fill=tk.BOTH, expand=True)

    welcome_tab = tk.Frame(notebook)
    notebook.add(welcome_tab, text="欢迎页")
    try:
        image_path = os.path.join(os.path.dirname(__file__), "hello.png")
        image = Image.open(image_path)
        image = image.resize((800, 600), Image.Resampling.LANCZOS)
        photo = ImageTk.PhotoImage(image)
        label = tk.Label(welcome_tab, image=photo)
        label.image = photo
        label.pack(pady=20)

        link = tk.Label(welcome_tab, text="访问 GitHub 项目主页", fg="blue", cursor="hand2", font=("Arial", 12, "underline"))
        link.pack(pady=10)
        link.bind("<Button-1>", lambda e: webbrowser.open(Github))

    except Exception as e:
        tk.Label(welcome_tab, text=f"无法加载欢迎图片：{e}").pack(pady=20)

    gwas_tab = gwas.create_frame(notebook)
    snp_tab = snpann.create_frame(notebook)
    primer_tab = primer.create_frame(notebook)

    notebook.add(gwas_tab, text="GWAS分析")
    notebook.add(snp_tab, text="SNP注释")
    notebook.add(primer_tab, text="引物/探针设计")

    root.mainloop()

if __name__ == "__main__":
    main()
