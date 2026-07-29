# DeepCaller

<p align="left"> 
  <img src="https://img.shields.io/badge/版本-1.0.0-blue" alt="版本">
  <img src="https://img.shields.io/badge/许可证-MIT-green" alt="许可证">
  <img src="https://img.shields.io/badge/python-3.9-blue" alt="Python">
  <img src="https://img.shields.io/badge/平台-linux-lightgrey" alt="平台">
</p>

**DeepCaller** 是一款基于深度学习的变异检测工具，专为多倍体基因组短读长数据中 SNP 和小片段 Indel 的精准检测而设计。它提供了五个针对四倍体和六倍体作物的预训练模型，并支持速度优先和性能优先两种推理模式。英文教程请参见 [English README](../README.md)。

> **注意**：本仓库配套论文目前正在审稿中，软件的完整使用权将在论文正式发表后开放。详情请参见 [LICENSE](../LICENSE)。

---

## 🏛️ 背景

<p align="center">
  <img src="flow.png" alt="DeepCaller 工作流程" width="800">
</p>

DeepCaller 的工作流程包含四个顺序步骤。**步骤一：** 对输入 BAM 文件进行过滤后，DeepCaller 逐位点分析比对数据，基于最小等位基因频率和测序深度的双重阈值筛选候选变异位点。**步骤二：** 将每个候选位点双链及其侧翼碱基编码为结构化的堆积张量（pileup tensor）。**步骤三：** 将张量输入由两层双向 LSTM（Bi-LSTM）和三层 ReLU 激活全连接层组成的循环神经网络（RNN），按倍性特定类别预测基因型（四倍体五类，六倍体七类）。**步骤四：** 根据预测基因型和比对数据生成 VCF 文件。

---

## 🌿 支持物种

| `--species`   | 常用名称           | 倍性   | 训练数据集       | 默认 |
|---------------|--------------------|--------|------------------|------|
| `potato`      | 四倍体马铃薯       | 四倍体 | C88              | ✓（四倍体） |
| `alfalfa`     | 苜蓿               | 四倍体 | Bolivia          | |
| `rose`        | 现代月季           | 四倍体 | Samantha         | |
| `sweetpotato` | 甘薯               | 六倍体 | Tanzania         | ✓（六倍体） |
| `syn_potato`  | 合成六倍体马铃薯   | 六倍体 | SyntheticPotato  | |

> 建议用户选择与目标物种最相近的物种模型；若不确定，推荐使用默认模型（四倍体默认使用 `potato`，六倍体默认使用 `sweetpotato`）。

---

## 🛠️ 安装

### 环境要求

- Linux (x86_64)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) ≥ 4.10

### 安装步骤

```bash
# 1. 克隆仓库
git clone https://github.com/JiaoLab2021/DeepCaller.git
cd DeepCaller

# 2. 创建并激活 conda 环境
conda env create -f deepcaller.yml
conda activate deepcaller

# 3. 安装 DeepCaller
pip install -e .

# 4. 验证安装
deepcaller --version
```

---

## 🚀 快速开始

`demo/` 目录中提供了一个小型演示数据集（第 10 号染色体 1 Mb 区域；四倍体马铃薯 C88）。

```bash
cd demo

deepcaller \
    -r DM8.1_chr10_100000_1100000.fa \
    -b C88_20x_chr10_100000_1100000.bam \
    -p 4 \
    --mode speed \
    -o demo_output.vcf
```

---

## 📖 使用说明

```
deepcaller -r <REF> -b <BAM> -p <PLOIDY> [options]
```

### 必需参数

| 参数 | 说明 |
|------|------|
| `-r`, `--ref` | 参考基因组 FASTA 文件 |
| `-b`, `--bam` | 输入 BAM 文件 |
| `-p`, `--ploidy` | 倍性水平：`4` 或 `6` |

### 输入/输出配置

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-o`, `--output` | `output.vcf` | 输出 VCF 文件（将进行 bgzip 压缩） |
| `-c`, `--chroms` | 全部 | 指定处理的染色体 |
| `--bed` | — | BED 文件，将变异检测限定于目标区域；设置后覆盖 `--chroms` |
| `-S`, `--sample` | `SAMPLE` | 显示在 VCF `#CHROM` 列头行中的样本名/ID |

### 处理选项

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-s`, `--species` | auto | 物种模型（参见"支持物种"） |
| `-m`, `--mode` | `speed` | 推理模式：`speed`（速度优先）或 `performance`（性能优先） |
| `-t`, `--cpus` | `24` | CPU 线程数；`-1` 表示使用全部可用线程 |
| `-d`, `--downsample` | 关闭 | 若测序深度超过目标深度（四倍体：50X，六倍体：80X），将 BAM 下采样至目标深度 |
| `--seed` | `42` | 用于 `samtools view -s` 下采样的随机种子 |
| `--min_af` | `0.10` | 候选位点最小等位基因频率 |
| `--rd_floor` | `10` | 候选位点最小测序深度 |

### 示例命令

```bash
# 四倍体马铃薯，全基因组，性能优先模式
deepcaller -r ref.fa -b sample.bam -p 4 --mode performance -o out.vcf -t 24

# 六倍体甘薯，指定染色体
deepcaller -r ref.fa -b sample.bam -p 6 -c chr1 chr2 chr3 -o out.vcf

# 苜蓿，仅分析目标区域（BED 文件）
deepcaller -r ref.fa -b sample.bam -p 4 --species alfalfa --bed targets.bed -o out.vcf

# 自定义样本名，对高深度染色体先下采样再检测变异
deepcaller -r ref.fa -b sample.bam -p 4 -S MySample -d -o out.vcf
```

---

## 📄 输出结果

DeepCaller 输出经 bgzip 压缩并建立 tabix 索引的 VCF 文件（`<output>.gz` 和 `<output>.gz.tbi`）。

### FORMAT 字段说明

| 字段 | 说明 |
|------|------|
| `GT` | 多倍体基因型（如四倍体单体型 `0/0/0/1`） |
| `GQ` | 基因型质量值 |
| `DP` | 该位点测序深度 |
| `AD` | 各等位基因深度（参考基因/替代基因） |
| `AF` | 等位基因频率 |

---

## 📝 引用

如果您在研究中使用了 DeepCaller，请引用：

> 

---

## ⚖️ 许可证

本项目采用 MIT 许可证，详情请参见 [LICENSE](../LICENSE)。  
软件的完整使用权将在配套论文正式发表后开放。

---

## 📬 联系方式

Kang Xiao · [xiaokangneuq@163.com](mailto:xiaokangneuq@163.com)
