# DNA-Sequencing-Result-Analyzer
A Python tool for analyzing DNA sequencing files (.ab1 and .fasta),  judging heterozygous/homozygous status and mutation type, and exporting results with parameters to CSV.

本工具用于 **自动分析 DNA 测序结果**，支持 **AB1 (测序原始文件)** 和 **FASTA (序列比对文件)**。  
能够帮助快速判断样本是否为 **杂合子 / 野生纯合 / 突变纯合**，并将分析结果和使用的参数 **导出到 CSV 文件**。

**小孩子不懂事，没事写着玩的.jpg**

---

## ✨ 功能特性
- 📂 支持输入：
  - `.dna` 原始 DNA 参考文件（野生型）
  - `.ab1` 测序峰图文件
  - `.fasta` 测序比对结果文件
- ⚙️ 一次性图形界面填写参数（阈值、裁剪长度等）
- 🔎 自动计算：
  - **AB1 文件** → 总碱基数、低质量比例、杂/纯合判定
  - **FASTA 文件** → 序列长度、匹配比例、野生/突变判定
- 🧬 综合判定规则：
  - AB1 判定为杂合 → **杂合子**
  - AB1 判定为纯合 → 根据 FASTA 判定：
    - **WT → 野生纯合**
    - **Mutant → 突变纯合**
- 📊 结果导出为 CSV：
  - 上方为各样本结果，可直接排序筛选
  - 下方记录运行参数，便于追溯

---

## ⚙️ 工作原理

1. **选择原始 DNA 文件 (.dna)**  
   - 作为参考序列（野生型）  
   - 自动兼容 UTF-8 / Latin1 编码  

2. **选择测序文件 (.ab1 + .fasta)**  
   - `.ab1`：提取碱基质量信息 (Phred 分数)，根据低质量碱基比例判断杂/纯合  
   - `.fasta`：将序列与野生型进行局部比对，计算匹配比例，判断 WT 或 Mutant  

3. **参数设置**  
   - 在一个弹窗中一次性填写：
     - 杂/纯合判断：去掉前多少 bp、质量阈值、低质量比例阈值  
     - 野生/突变判断：两端裁剪长度、匹配比例阈值  

4. **结果判定逻辑**  
   - AB1 判定为 "可能杂合" → 样本为 **杂合子**  
   - 否则根据 FASTA 判定：
     - 匹配比例 ≥ 阈值 → **野生纯合**  
     - 否则 → **突变纯合**  

5. **结果导出 CSV**  
   - 每个样本一行，包含：
     - `AB1_总碱基数`、`AB1_低质量比例`、`FASTA_匹配比例`、综合判定  
   - 文件最后附带本次运行所使用的参数  

---

## 📥 安装与运行

### 依赖
- Python 3.8+
- Biopython
- Tkinter（Python 自带）

### 安装依赖
```bash
pip install biopython
