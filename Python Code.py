from tkinter import Tk, filedialog, simpledialog
import tkinter as tk
from Bio import SeqIO, Align
import os
import csv

def get_params():
    """
    弹出一个窗口让用户一次性填写所有参数
    返回字典 params
    """

    def confirm():
        nonlocal params
        params = {
            "start_pos": int(entry_start_pos.get() or 30),
            "q_threshold": int(entry_q_threshold.get() or 20),
            "low_quality_ratio_threshold": float(entry_low_quality_ratio_threshold.get() or 0.15),
            "match_threshold": float(entry_match_threshold.get() or 0.85),
            "trim_bp": int(entry_trim_bp.get() or 15),
        }
        win.destroy()

    params = None
    win = tk.Toplevel()
    win.title("参数设置")

    # 字体样式
    hint_font = ("Arial", 8, "italic")

    # --- AB1 参数 ---
    tk.Label(win, text="杂/纯合判断中 去掉前多少 bp:").grid(row=0, column=0, sticky="e")
    entry_start_pos = tk.Entry(win)
    entry_start_pos.insert(0, "30")
    entry_start_pos.grid(row=0, column=1)
    tk.Label(win, text="建议去掉前端低质量片段 (默认 30)", font=hint_font, fg="gray").grid(row=1, column=0, columnspan=2, sticky="w")

    tk.Label(win, text="杂/纯合判断中 Phred 质量阈值:").grid(row=2, column=0, sticky="e")
    entry_q_threshold = tk.Entry(win)
    entry_q_threshold.insert(0, "20")
    entry_q_threshold.grid(row=2, column=1)
    tk.Label(win, text="一般取 20 ，正确率较高 (默认 20)", font=hint_font, fg="gray").grid(row=3, column=0, columnspan=2, sticky="w")

    tk.Label(win, text="杂/纯合判断中 低质量比例阈值:").grid(row=4, column=0, sticky="e")
    entry_low_quality_ratio_threshold = tk.Entry(win)
    entry_low_quality_ratio_threshold.insert(0, "0.15")
    entry_low_quality_ratio_threshold.grid(row=4, column=1)
    tk.Label(win, text="低质量碱基比例超过阈值则判定为杂合 (默认 0.15)", font=hint_font, fg="gray").grid(row=5, column=0, columnspan=2, sticky="w")

    # --- FASTA 参数 ---
    tk.Label(win, text="野生/突变判断中 匹配比例阈值:").grid(row=6, column=0, sticky="e")
    entry_match_threshold = tk.Entry(win)
    entry_match_threshold.insert(0, "0.85")
    entry_match_threshold.grid(row=6, column=1)
    tk.Label(win, text="比对匹配比例 ≥ 阈值判定为野生型 (默认 0.85)", font=hint_font, fg="gray").grid(row=7, column=0, columnspan=2, sticky="w")

    tk.Label(win, text="野生/突变判断中 去掉首尾 bp 数:").grid(row=8, column=0, sticky="e")
    entry_trim_bp = tk.Entry(win)
    entry_trim_bp.insert(0, "15")
    entry_trim_bp.grid(row=8, column=1)
    tk.Label(win, text="去掉测序两端的低质量片段 (默认 15)", font=hint_font, fg="gray").grid(row=9, column=0, columnspan=2, sticky="w")

    # --- 确认按钮 ---
    tk.Button(win, text="确定", command=confirm, bg="lightblue").grid(row=10, column=0, columnspan=2, pady=10)

    win.grab_set()  # 阻塞主线程
    win.wait_window()

    return params

# ----------------------------
# 1. 选择原始 DNA 文件 (.dna)
# ----------------------------
Tk().withdraw()
wt_file_path = filedialog.askopenfilename(filetypes=[("DNA files", "*.dna")])
if not wt_file_path:
    print("未选择原始 DNA 文件，程序退出")
    exit()

# 尝试不同编码读取
try:
    with open(wt_file_path, "r", encoding="utf-8") as f:
        wt_sequence = f.read().replace("\n", "").upper()
except UnicodeDecodeError:
    with open(wt_file_path, "r", encoding="latin1", errors="ignore") as f:
        wt_sequence = f.read().replace("\n", "").upper()

print("选择的原始 DNA 文件:", wt_file_path)

# ----------------------------
# 2. 选择测序文件 (.ab1 + .fasta)
# ----------------------------
seq_file_paths = filedialog.askopenfilenames(filetypes=[("Sequencing files", "*.ab1 *.fasta")])
if not seq_file_paths:
    print("未选择测序文件，程序退出")
    exit()

print(f"选择的测序文件数量: {len(seq_file_paths)}")

# ----------------------------
# 3. 参数设置（一个窗口输入）
# ----------------------------
params = get_params()

start_pos = params["start_pos"]
q_threshold = params["q_threshold"]
low_quality_ratio_threshold = params["low_quality_ratio_threshold"]
match_threshold = params["match_threshold"]
trim_bp = params["trim_bp"]

print(f"使用参数：AB1 去掉前 {start_pos}, Q阈值={q_threshold}, 低质量比例阈值={low_quality_ratio_threshold}")
print(f"使用参数：FASTA 匹配阈值={match_threshold}, 去掉首尾 {trim_bp}bp")

# ----------------------------
# 4. 处理 AB1 文件
# ----------------------------
ab1_results = {}

for file_path in seq_file_paths:
    if not file_path.endswith(".ab1"):
        continue

    record = SeqIO.read(file_path, "abi")
    quality = record.letter_annotations["phred_quality"]

    # 去掉前 N bp
    qual_trimmed = quality[start_pos:]
    total_bases = len(qual_trimmed)
    low_quality_count = sum(1 for q in qual_trimmed if q < q_threshold)
    low_quality_ratio = low_quality_count / total_bases if total_bases > 0 else 0

    if low_quality_ratio >= low_quality_ratio_threshold:
        judgement = "可能杂合"
    else:
        judgement = "可能纯合"

    sample_name = os.path.splitext(os.path.basename(file_path))[0]
    ab1_results[sample_name] = {
        "AB1_总碱基数": total_bases,
        "AB1_低质量碱基数": low_quality_count,
        "AB1_低质量比例": f"{low_quality_ratio:.2f}",
        "AB1_判断": judgement,
    }

# ----------------------------
# 5. 处理 FASTA 文件
# ----------------------------
aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1
aligner.mismatch_score = -100
aligner.open_gap_score = -10
aligner.extend_gap_score = -10

fasta_results = {}

for file_path in seq_file_paths:
    if not file_path.endswith(".fasta"):
        continue

    record = SeqIO.read(file_path, "fasta")
    seq = str(record.seq).upper()

    # 去掉首尾
    if len(seq) > 2 * trim_bp:
        seq = seq[trim_bp:-trim_bp]
    else:
        seq = ""

    if len(seq) > 0:
        score = aligner.score(wt_sequence, seq)
        match_ratio = score / len(seq)
    else:
        score = 0
        match_ratio = 0

    judgement = "WT" if match_ratio >= match_threshold else "Mutant"

    sample_name = os.path.splitext(os.path.basename(file_path))[0]
    fasta_results[sample_name] = {
        "FASTA_序列长度": len(seq),
        "FASTA_匹配比例": f"{match_ratio:.2f}",
        "FASTA_判断": judgement,
    }

# ----------------------------
# 6. 合并结果 & 综合判定
# ----------------------------
results = []
all_samples = set(list(ab1_results.keys()) + list(fasta_results.keys()))

for sample in all_samples:
    ab1_res = ab1_results.get(sample, {})
    fasta_res = fasta_results.get(sample, {})

    # 综合判定规则：
    # 如果 AB1 判定为杂合 → 杂合子
    # 否则根据 FASTA 判断：WT=野生纯合，Mutant=突变纯合
    if ab1_res.get("AB1_判断") == "可能杂合":
        final_judgement = "杂合子"
    else:
        if fasta_res.get("FASTA_判断") == "WT":
            final_judgement = "野生纯合"
        elif fasta_res.get("FASTA_判断") == "Mutant":
            final_judgement = "突变纯合"
        else:
            final_judgement = "未知"

    merged = {"样本": sample}
    merged.update(ab1_res)
    merged.update(fasta_res)
    merged["综合判定"] = final_judgement

    results.append(merged)

# ----------------------------
# 7. 保存结果 CSV（结果在前，参数在后）
# ----------------------------

# 先按照样本排序
results = sorted(results, key=lambda x: x["样本"])

output_folder = os.path.dirname(seq_file_paths[0])
output_file = os.path.join(output_folder, "sequencing_combined_results.csv")

# 需要记录的参数
params = {
    "原始DNA文件": os.path.basename(wt_file_path),
    "杂/纯合判断中 起始去掉碱基数": start_pos,
    "杂/纯合判断中 低质量碱基判定阈值": q_threshold,
    "杂/纯合判断中 低质量比例阈值": low_quality_ratio_threshold,
    "野生/突变判断中 两端裁剪碱基数": trim_bp,
    "野生/突变判断中 FASTA匹配比例阈值": match_threshold,
}

with open(output_file, mode="w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)

    # 写结果表头
    fieldnames = ["样本", 
                  "AB1_总碱基数", "AB1_低质量碱基数", "AB1_低质量比例", "AB1_判断",
                  "FASTA_序列长度", "FASTA_匹配比例", "FASTA_判断",
                  "综合判定"]
    writer.writerow(fieldnames)

    # 写结果数据
    for row in results:
        writer.writerow([row.get(fn, "") for fn in fieldnames])

    # 空行分隔
    writer.writerow([])
    writer.writerow(["参数设置"])
    for k, v in params.items():
        writer.writerow([k, v])

print(f"综合结果已保存到: {output_file}")
