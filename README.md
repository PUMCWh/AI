# CKM 论文分析代码整理版

本仓库已按“可交付流程”拆分为 3 个主阶段，并将探索性代码单独隔离。

## 目录结构

- `01_data_processing/`：数据整理与处理
- `02_analysis/`：数据分析（含初步图表）
- `03_publication_outputs/`：结果展示（出版级表格和图片）生成
- `99_exploratory/`：探索性/参考性脚本（不纳入主交付流水线）

---

## 推荐执行顺序

### 1) 数据整理与处理（`01_data_processing/`）
建议顺序：
1. `数据准备.R`
2. `偏相关网络构建.R`
3. `社区发现.R`
4. `社区映射回个体.R`
5. `主导社区分配.R`
6. `QC_Step4_快速检查.R`（质控）

### 2) 数据分析（含初步图表，`02_analysis/`）
建议顺序：
1. `cox回归分析.R`
2. `亚组与交互作用分析.R`
3. `限制性立方样条分析.R`
4. `人群归因分析.R`
5. `ICD-10 Chapter Enrichment Analysis for CKM Disease Communities.R`
6. `Node Role Classification (Guimerà-Amaral framework).R`
7. `Normalized Connectivity Index (NCI) Matrix across 16 Disease Communities.R`
8. `Bootstrap C-statistic and ΔC 95% Confidence Intervals.R`
9. `Top1-share distribution with τ = 0.40 threshold justification.R`
10. `Leiden Resolution Sweep.R`
11. `社区ICD章节富集分析.R`
12. `社区间标准化连接强度热图.R`
13. `桥接疾病的节点角色分类.R`

### 3) 结果展示（出版级表格和图片，`03_publication_outputs/`）
建议顺序：
1. `生成正文表1.R`
2. `生成正文表2.R`
3. `生成正文图3.R`
4. `生成出版级森林图（基于cox回归）.R`
5. `STROBE Flow Diagram for Figure.R`
6. `Fig S3 Top 100 diseases by Multimorbidity Centrality (MMC), circular bar chart.R`
7. `出版级分辨率扫描结果图.R`
8. `Fig 3 Kaplan-Meier 生存曲线 (Publication grade BMC Medicine 投稿版.R`
9. `Generate Supplementary Tables S2 - S9 in publication format.R`
10. `# Step 15-16-17 三连发 — 补表 S8、补表 S9、补图 S7.R`

---

## 说明

- 本次整理优先做“交付结构重构”，未改动各脚本内部统计逻辑。
- 若你愿意，我下一步可以继续：
  1) 统一输入/输出路径配置（避免脚本内硬编码路径）；
  2) 做一个一键主控脚本（按阶段自动执行）；
  3) 标注每个脚本的“必需输入文件”和“产出文件清单”。


## 可交付一键脚本（新）

为避免改动原分析与绘图参数，本仓库新增 `deliverable_pipeline/`：

- `deliverable_pipeline/01_data_processing_pipeline.R`
- `deliverable_pipeline/02_analysis_pipeline.R`
- `deliverable_pipeline/03_publication_outputs_pipeline.R`

它们按阶段顺序直接调用原脚本（`source(..., chdir=TRUE)`），因此图形参数与统计逻辑保持不变。

建议执行：

```bash
Rscript deliverable_pipeline/01_data_processing_pipeline.R
Rscript deliverable_pipeline/02_analysis_pipeline.R
Rscript deliverable_pipeline/03_publication_outputs_pipeline.R
```

## 全新重写版代码（不调用旧脚本）

如果你希望完全不依赖历史脚本、直接使用“新生成的可交付分析/出图代码”，请使用：

- `regenerated/02_analysis_and_prelim_plots_regenerated.R`
- `regenerated/03_publication_tables_figures_regenerated.R`

特点：
- 不 `source()` 旧脚本；
- 直接读取 `CKM_COX_DATA` 指向的数据文件；
- 生成新的分析结果表和图。

示例：

```bash
export CKM_COX_DATA=data/cox_dataset_with_dominant_community.csv
Rscript regenerated/02_analysis_and_prelim_plots_regenerated.R
Rscript regenerated/03_publication_tables_figures_regenerated.R
```
