# DeepCaller

<p align="left"> 
  <img src="https://img.shields.io/badge/版本-1.0.1-blue" alt="版本">
  <img src="https://img.shields.io/badge/许可证-MIT-green" alt="许可证">
  <img src="https://img.shields.io/badge/python-3.9-blue" alt="Python">
  <img src="https://img.shields.io/badge/平台-linux-lightgrey" alt="平台">
</p>

**DeepCaller** 是一款基于深度学习的变异检测工具，用于从短读长数据中精准检测多倍体基因组的 SNP 和小片段 Indel。它提供五个针对四倍体和六倍体作物的预训练模型，推理阶段可在普通 CPU 上运行。英文教程请参见 [English README](../README.md)。

> **注意**：本仓库配套论文正在审稿中，软件的完整使用权将在论文正式发表后开放。详情请参见 [LICENSE](../LICENSE)。

---

## 🏛️ 背景

<p align="center">
  <img src="flow.png" alt="DeepCaller 工作流程" width="800">
</p>

DeepCaller 的工作流程包含四个顺序步骤。**步骤一，候选发现：** 对输入 BAM 文件进行过滤后，逐位点扫描，基于替代等位基因频率与测序深度的双重阈值筛选候选变异位点。**步骤二，读段分组：** 将覆盖每个候选位点的读段按其支持的替代等位基因分组，分组数不超过样本倍性。**步骤三，特征编码：** 将每个等位基因分组的堆积（pileup）连同其侧翼位点编码为形状为 (2*w* + 1) × 15 的结构化张量。**步骤四，剂量预测：** 由权重共享的 LSTM 汇总每个分组，跨组自注意力在各组间交换上下文信息，解码器在硬性倍性预算约束下自回归地预测每个候选等位基因的拷贝数，并据此生成 VCF 文件。

---

<a id="支持物种"></a>
## 🌿 支持物种

| `--species`             | 常用名称           | 倍性   | 训练数据集       | 默认         |
|-------------------------|--------------------|--------|------------------|--------------|
| `C88_Potato`            | 四倍体马铃薯       | 四倍体 | C88              | ✓（四倍体）  |
| `Bolivia_Alfalfa`       | 苜蓿               | 四倍体 | Bolivia          | |
| `Samantha_Rose`         | 现代月季           | 四倍体 | Samantha         | |
| `SyntheticPotato_Potato`| 合成六倍体马铃薯   | 六倍体 | 合成六倍体       | ✓（六倍体）  |
| `Tanzania_Sweetpotato`  | 甘薯               | 六倍体 | Tanzania         | |

> 若不指定 `--species`，DeepCaller 将按倍性使用默认模型（四倍体为 `C88_Potato`，六倍体为 `SyntheticPotato_Potato`）。默认模型跨物种泛化良好，对没有专用模型的基因组是稳妥选择；当在意几个百分点的精度差异时，建议在一小段区域上比较各候选模型后再定。

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
    -o demo_output.vcf
```

---

```bash
## 📖 使用说明
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
| `-o`, `--out` | `output.vcf` | 输出 VCF 文件（bgzip 压缩并建立 tabix 索引） |
| `-c`, `--chroms` | 全部 | 指定处理的染色体 |
| `-l`, `--bed` | — | BED 文件，将变异检测限定于目标区域；设置后忽略 `--chroms` |
| `--sample` | `SAMPLE` | 写入 VCF `#CHROM` 表头行的样本名/ID |
| `--work_dir` | 自动 | 临时目录；默认在工作目录下新建唯一目录，成功后自动删除 |
| `--keep_tmp` | 关闭 | 保留临时目录以便调试 |

### 处理选项

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--species` | 随倍性而定 | 物种模型（参见[支持物种](#支持物种)）；`-p 4` 默认 `C88_Potato`，`-p 6` 默认 `SyntheticPotato_Potato` |
| `-t`, `--cpus` | `24` | CPU 线程数；`-1` 表示使用全部可用线程 |
| `--downsample` | 关闭 | 当全基因组深度超过目标深度（四倍体 50×，六倍体 80×）时，将 BAM 下采样至目标深度 |
| `--seed` | `42` | 用于 `samtools view -s` 下采样的随机种子 |
| `-v`, `--min_af` | `0.1` | 候选等位基因的最小替代等位基因频率 |
| `-d`, `--rd_floor` | `8` | 候选位点的最小测序深度 |
| `--min_mq` | `5` | 堆积时保留读段的最小比对质量 |
| `--max_id_len` | `50` | 作为候选等位基因的最大 indel 长度 |
| `--batch_size` | `8192` | 模型推理批大小 |

### 示例命令

```bash
# 四倍体马铃薯，全基因组，24 线程
deepcaller -r ref.fa -b sample.bam -p 4 -o out.vcf -t 24

# 六倍体甘薯模型，指定染色体
deepcaller -r ref.fa -b sample.bam -p 6 --species Tanzania_Sweetpotato -c chr1 chr2 chr3 -o out.vcf

# 苜蓿，仅分析目标区域（BED 文件）
deepcaller -r ref.fa -b sample.bam -p 4 --species Bolivia_Alfalfa -l targets.bed -o out.vcf

# 自定义样本名；对高深度数据先下采样再检测变异
deepcaller -r ref.fa -b sample.bam -p 4 --sample MySample --downsample -o out.vcf
```

---

## 📄 输出结果

DeepCaller 输出经 bgzip 压缩并建立 tabix 索引的 VCF 文件（`<output>.gz` 和 `<output>.gz.tbi`）。

### FORMAT 字段说明

| 字段 | 说明 |
|------|------|
| `GT` | 多倍体基因型（如四倍体单拷贝位点 `0/0/0/1`） |
| `GQ` | 基因型质量值 |
| `DP` | 该位点测序深度 |
| `AD` | 各等位基因深度（参考，替代） |
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
