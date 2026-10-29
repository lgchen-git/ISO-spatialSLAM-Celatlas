# 长读长空间转录组与 SLAM 分析流程

## 1. 流程目的

本流程用于处理带空间条形码的长读长转录组数据，完成原始数据质控、空间条形码识别、剪接比对、全长转录本发现、SLAM 分子分类、长读长空间表达聚类以及空间表达与 SLAM 联合分析。

流程由配置文件和样本表驱动，统一入口为 `workflow/run.sh`。所有科学参数均由用户在配置文件中明确填写，运行过程中不执行参数网格扫描或自动参数推荐，也不需要交互式或 AI 辅助。

## 2. 输入文件

每个样本至少需要以下输入：

| 输入 | 说明 |
|---|---|
| FASTQ | 长读长原始测序数据 |
| 空间条形码白名单 | 芯片有效条形码及空间坐标 |
| 参考基因组 FASTA | 用于剪接比对和碱基替换统计 |
| 基因注释 GTF | 用于基因、外显子和转录本注释 |
| 样本表 | 记录样本 ID、芯片 ID、物种、化学版本和输入路径 |
| 项目配置 | 记录软件路径、过滤阈值、聚类参数和报告参数 |

可选输入包括预构建的 minimap2 索引、剪接位点 BED、外部变异 VCF 和已有的 FLAMES/gffcompare 结果目录。

## 3. 八阶段分析

| 阶段 | 分析内容 | 主要软件或方法 | 核心产物 |
|---|---|---|---|
| `01_read_qc` | 统计读段数、碱基数、读长、N50 和质量值 | NanoPlot、流程内汇总脚本 | 质控表、读长图、质量图 |
| `02_barcode_cdna` | 识别接头和 linker，匹配空间条形码，提取 cDNA | 流程内结构扫描与白名单匹配 | 条形码状态、cDNA FASTQ、空间坐标 |
| `03_align` | 将 cDNA 剪接比对到参考基因组，并执行质量门控 | minimap2、samtools | 排序 BAM、索引、比对统计 |
| `04_full_length_transcriptome` | 发现候选全长转录本并与参考注释比较 | FLAMES、gffcompare | GFF3、FASTA、类别表、结构图 |
| `05_slam` | 汇总分子层面的 T→C/A→G 等替换证据，按手工阈值分类 | 流程内分子调用与统计 | SLAM 分子表、基因汇总、12 类碱基替换图 |
| `06_whole_transcriptome` | 构建 Bin×gene 长读长空间表达矩阵并聚类 | Scanpy、Leiden、空间邻接图 | H5AD、空间簇、UMAP、Marker、空间连续性 |
| `07_slam_joint` | 将 SLAM 分子映射到固定空间 Bin 和表达簇 | 流程内联合统计 | 簇级和基因×簇级 SLAM 表及图 |
| `08_report` | 汇总所有阶段，生成中文离线报告 | 流程内 HTML 报告生成器 | `<sample_id>_report.html`、交付清单 |

内部阶段名 `06_whole_transcriptome` 为稳定接口；报告中称为“长读长空间转录组”，表示由全部有效长读长分子构建的空间基因表达矩阵。

`06_whole_transcriptome` 另外输出 bin50 空间 QC 表，包括每个 Bin 的 UMI、长读长分子数、detected gene 数和探索性 dynamic proxy；样本级 isoform 数来自 FLAMES，不代表每个 Bin 的 isoform 数。

## 4. 关键分析方法

### 4.1 全长转录组

FLAMES 从长读长剪接比对结果中发现候选异构体，gffcompare 将候选转录本与参考注释进行结构比较。输出包括参考一致转录本、新剪接组合、包含关系和其他结构类别。

该阶段侧重转录本结构。若 gffcompare 的 TMAP 中表达字段为零，不将其解释为可靠的 TPM/FPKM 定量。

### 4.2 长读长空间表达聚类

空间表达阶段按配置的 Bin 大小聚合分子，过滤低覆盖 Bin 和低检出基因，然后执行：

1. 总量归一化与 `log1p` 转换；
2. 高变基因选择和标准化；
3. PCA 与测序深度线性校正；
4. 构建表达邻接图；
5. 以配置权重加入空间四邻域图；
6. 使用手工指定分辨率进行 Leiden 聚类；
7. 合并小于配置阈值的簇；
8. 生成 UMAP、空间簇图、Marker 和空间连续性统计。

old/new 分子不分别聚类，SLAM 类别不会参与表达簇的定义。

### 4.3 SLAM 联合分析

SLAM 分子按空间坐标落入与表达聚类一致的 Bin，然后在以下层级汇总：

- 空间 Bin；
- 表达簇；
- 基因；
- 基因×表达簇。

新生支持比例以 `new_supported + no_conversion_detected` 为信息分子分母，并报告 Wilson 95% 置信区间。当前分子表不含 `read→isoform` 映射，因此不输出异构体级 SLAM 结果。

正常分析在 `08_report` 结束。`09_promo_figures` 是报告后的可选交付模块，使用 `--promo-figures` 显式启动，可生成单面板 PPT 图、候选基因案例和 read-to-isoform 空间探索。

## 5. 配置与运行

准备项目配置和样本表：

```bash
cp conf/lr_celatlas_spatial_slam.template.yaml conf/<project>.yaml
cp conf/samples.lr_celatlas_spatial_slam.template.tsv conf/<project>.samples.tsv
```

填写所有 `CHANGE_ME` 项及科学参数后，执行：

```bash
make validate CONFIG=conf/<project>.yaml
bash workflow/run.sh --config conf/<project>.yaml --dry-run
bash workflow/run.sh --config conf/<project>.yaml
```

从指定阶段继续运行：

```bash
bash workflow/run.sh --config conf/<project>.yaml \
  --from 05_slam --to 08_report
```

只重新生成报告：

```bash
bash workflow/run.sh --config conf/<project>.yaml \
  --from 08_report --to 08_report --force-stage 08_report
```

仅使用既有科学结果重绘图形并更新报告：

```bash
bash workflow/run.sh --config conf/<project>.yaml --refresh-visuals
```

## 6. 结果结构与可追溯性

每个样本使用独立的八阶段结果目录。每个阶段下包含：

| 子目录 | 内容 |
|---|---|
| `01_tables` | 可交付统计表 |
| `02_figures` | PNG/PDF 图形 |
| `03_logs` | 命令和运行日志 |
| `04_tool_output` | 第三方软件原始输出 |
| `05_metadata` | 软件版本、Git 信息、输入输出清单和完成标记 |

阶段完成标记记录项目配置的 SHA-256。配置发生变化时，流程不会静默复用不匹配的完成标记。

最终报告位置为：

```text
<results_root>/<sample_id>/08_report/<sample_id>_report.html
```

## 7. 验证

代码或配置变更后至少执行：

```bash
make lint
make smoke
python3 scripts/validate_delivery_report.py \
  --report <results_root>/<sample_id>/08_report/<sample_id>_report.html \
  --manifest <results_root>/<sample_id>/08_report/05_metadata/delivery_manifest.tsv
```

验收内容包括脚本语法、配置与样本表契约、八阶段顺序、结果目录结构、报告章节、离线资源、表格分页和最终报告唯一性。
