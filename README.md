# 胃癌腹膜转移空间分析

公开转录组上的计算结果。GSE251950 的 9 张原发灶 Visium 没有留下成纤维–巨噬界面超额。GSE183904 里两条 C3 × C3AR1 样本乘积在 3 份腹膜肿瘤高于原发肿瘤，升高主要来自 C3。GSE308231 没有复现这个乘积。

## 文件

- 中文分析报告：`latex/analysis/main.pdf`
- 英文稿：`latex/mdpi/manuscript.docx` 与 `latex/mdpi/manuscript.pdf`
- 上皮状态轴是另一篇稿：`latex/manuscript/main.pdf`（GSE249874）
- 客户包：`胃癌腹膜转移分析_delivery/`，其中报告只有中文分析 PDF
- 结果表：`06_真实分析结果/`
- 脚本：`03_脚本/`

## 数据

原始计数矩阵不在本仓库。正式统计在分析服务器 `shengxin_01`：空间分析目录 `/root/autodl-tmp/gc_spatial_niche`，状态轴目录 `/root/autodl-tmp/gc_state_shift`。

公开来源：GSE251950、GSE183904、GSE308231、GSE163558、GSE228598、GSE249874。

克隆后不能重跑原始矩阵。出图脚本读取已保存的结果表。

## 限制

腹膜乘积是看过类均值之后指定的 7 条，q 只在这 7 条内部计算。3 份腹膜肿瘤不能写成已经验证的免疫抑制轴。没有湿实验、生存分析，也没有细胞级空间分辨率。
