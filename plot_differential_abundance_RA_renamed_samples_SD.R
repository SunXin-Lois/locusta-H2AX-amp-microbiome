## ============================================================
## 差异丰度图 + 相对丰度柱状图
## 数据：f_Enterobacteriaceae 内仅保留 Cronobacter、Enterobacter、
##       Trabulsiella；其余属及 Unclassified 一并并入 Unclassified
## RA 计算方式：属 reads / 全部样本总 reads（非仅肠杆菌科）
## ============================================================

.libPaths(c(.libPaths(), Sys.getenv("R_LIBS_USER")))

library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
# ── 用户参数 ────────────────────────────────────────────────────────
# ★ 修改此处指定输入文件路径
INPUT_FILE <- "table_filtered_w_tax_订正肠杆菌科分类_改分组_改样本名.tsv"

# ★ 仅单独展示的属；其余肠杆菌科属全部并入 Unclassified
KEEP_GENERA <- c("Cronobacter", "Enterobacter", "Trabulsiella")

# ★ 对照组，以及组名与对应的列名前缀（顺序决定图例顺序）
CONTROL_GROUP  <- "PTW"
GROUP_PREFIXES <- c("PTW", "1 dpi", "3 dpi", "5 dpi")

# ★ 颜色（与 GROUP_PREFIXES 对应）
GROUP_COLORS <- c("PTW" = "#B5B5B6", "1 dpi" = "#DFA213",
                  "3 dpi" = "#57ABDC", "5 dpi" = "#CC5D17")

# ★ 输出文件名（NULL = 只预览不保存）
OUTPUT_PDF <- "differential_abundance_RA_renamed_samples_mean_SD_3genera_260819.pdf"
OUTPUT_GENUS_COUNTS <- "Enterobacteriaceae_genus_merged_3genera_counts.tsv"
OUTPUT_GENUS_RA <- "Enterobacteriaceae_genus_merged_3genera_RA.tsv"
OUTPUT_RA_TABLE <- "differential_abundance_RA_renamed_samples_mean_SD_percent_by_timepoint_3genera.tsv"
OUTPUT_RA_TABLE_WIDE <- "differential_abundance_RA_renamed_samples_mean_SD_percent_wide_3genera.tsv"
OUTPUT_RA_TABLE_FORMATTED <- "differential_abundance_RA_renamed_samples_mean_SD_percent_formatted_3genera.tsv"
OUTPUT_WIDTH  <- 4.5   # 英寸
OUTPUT_HEIGHT <- 4.7   # 英寸

# ★ 右图误差棒粗细（linewidth）和帽宽（width）
ERRORBAR_LINEWIDTH <- 0.125    # ← 增大使误差棒更粗
ERRORBAR_CAPWIDTH  <- 0.25   # ← 增大使帽子更宽

# ★ 右图条形宽度与躲避宽度
DODGE_WIDTH <- 0.85
BAR_WIDTH   <- 0.70


# ── 读取数据 ────────────────────────────────────────────────────────
df_raw <- read.table(INPUT_FILE, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE, quote = "")

# 识别样本列（非分类列）
tax_cols    <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
sample_cols <- setdiff(colnames(df_raw), tax_cols)


# ── 计算相对丰度（RA = 属 reads / 所有样本总 reads）────────────────
# 先用全体数据（所有科）计算每样本总 reads，作为归一化分母
total_reads <- colSums(df_raw[, sample_cols])

# 筛选肠杆菌科
df_ent <- df_raw[df_raw$Family == "f_Enterobacteriaceae", ]

# 清理属名并合并：
#   仅保留 Cronobacter / Enterobacter / Trabulsiella
#   其余属（含 Unclassified 及 Klebsiella 等）全部并入 Unclassified
df_ent$Display_Genus <- sub("^g_", "", df_ent$Genus)
df_ent$Display_Genus[df_ent$Display_Genus == "Unclassified_f_Enterobacteriaceae"] <- "Unclassified"
df_ent$Display_Genus[df_ent$Display_Genus == "Pseudomonas_f_Enterobacteriaceae"] <- "Pseudomonas"
df_ent$Display_Genus[!df_ent$Display_Genus %in% KEEP_GENERA] <- "Unclassified"

# 按清理后的属名合并 reads
genus_counts <- aggregate(df_ent[, sample_cols],
                          by  = list(Genus = df_ent$Display_Genus),
                          FUN = sum)
rownames(genus_counts) <- genus_counts$Genus
genus_counts$Genus <- NULL

# 每个属在每个样本的 RA（相对于全体样本总 reads）
genus_ra <- sweep(genus_counts, 2, total_reads, "/")

# 展示顺序：三个指定属 + Unclassified（全部保留，不再做 Top N 筛选）
display_order <- c(KEEP_GENERA, "Unclassified")
present_genera <- intersect(display_order, rownames(genus_ra))
genus_ra_top <- genus_ra[present_genera, , drop = FALSE]

if (!is.null(OUTPUT_GENUS_COUNTS)) {
  genus_counts_out <- genus_counts[present_genera, , drop = FALSE]
  write.table(
    cbind(Genus = rownames(genus_counts_out), genus_counts_out),
    file = OUTPUT_GENUS_COUNTS,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("合并后属水平 counts 表已保存: ", OUTPUT_GENUS_COUNTS)
}

if (!is.null(OUTPUT_GENUS_RA)) {
  genus_ra_out <- genus_ra[present_genera, , drop = FALSE]
  write.table(
    cbind(Genus = rownames(genus_ra_out), genus_ra_out),
    file = OUTPUT_GENUS_RA,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("合并后属水平 RA 表已保存: ", OUTPUT_GENUS_RA)
}


# ── 统计：各组均值、SD、Log2FC vs 对照组 ──────────────────────────
# 识别每组的样本列
group_cols <- lapply(GROUP_PREFIXES, function(pfx) {
  grep(paste0("^", pfx, "-"), colnames(genus_ra_top), value = TRUE)
})
names(group_cols) <- GROUP_PREFIXES

if (any(lengths(group_cols) == 0)) {
  missing_groups <- names(group_cols)[lengths(group_cols) == 0]
  stop("未识别到以下分组的样本列: ", paste(missing_groups, collapse = ", "))
}

results_list <- list()

for (genus in rownames(genus_ra_top)) {
  ra_vals <- as.numeric(genus_ra_top[genus, group_cols[[CONTROL_GROUP]]])
  m_ck    <- mean(ra_vals)
  sd_ck   <- sd(ra_vals)

  # 对照组行（L2FC = 0）
  results_list[[length(results_list) + 1]] <- data.frame(
    Taxon = genus, Group = CONTROL_GROUP,
    Mean = m_ck, SD = sd_ck, L2FC = 0
  )

  # 各处理组
  for (gn in setdiff(GROUP_PREFIXES, CONTROL_GROUP)) {
    vals  <- as.numeric(genus_ra_top[genus, group_cols[[gn]]])
    m_tr  <- mean(vals)
    sd_tr <- sd(vals)
    l2fc  <- log2((m_tr + 1e-6) / (m_ck + 1e-6))

    results_list[[length(results_list) + 1]] <- data.frame(
      Taxon = genus, Group = gn,
      Mean = m_tr, SD = sd_tr, L2FC = l2fc
    )
  }
}

plot_data <- do.call(rbind, results_list)

# Y 轴顺序：按对照组均值升序（最大值在最上方）
taxon_order <- plot_data %>%
  filter(Group == CONTROL_GROUP) %>%
  arrange(Mean) %>%            # 升序 → coord_flip 后顶部最大
  pull(Taxon)

plot_data$Taxon <- factor(plot_data$Taxon, levels = taxon_order)
plot_data$Group <- factor(plot_data$Group, levels = GROUP_PREFIXES)

# ── 导出统计表：各时间点 RA 百分比（Top N 属）────────────────────────
group_n <- sapply(group_cols, length)

ra_percent_table <- plot_data %>%
  mutate(
    Taxon = as.character(Taxon),
    Group = as.character(Group),
    N = as.integer(group_n[Group]),
    RA_percent_mean = Mean * 100,
    RA_percent_sd = SD * 100
  ) %>%
  select(Taxon, Group, N, RA_percent_mean, RA_percent_sd) %>%
  arrange(factor(Taxon, levels = taxon_order), factor(Group, levels = GROUP_PREFIXES))


# ── P1：左图 — Log2 Fold Change 点图 ───────────────────────────────
# 使用 shape = 21（可同时设置填充色 fill 和边线色 color 的实心圆）
# aes 中用 fill = Group 控制各组填充色，color 固定为黑色边线
# stroke 控制黑色边线的粗细  ← 修改 stroke 调整边线粗细
p1 <- ggplot(
    subset(plot_data, Group != CONTROL_GROUP),
    aes(x = L2FC, y = Taxon, fill = Group)   # 改为 fill 映射分组颜色
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#B5B5B6") +
  geom_point(
    shape  = 21,      # 带边线的实心圆；22=方形 23=菱形 24=三角，均支持 fill+color
    size   = 3,     # 点的大小          ← 修改此处
    stroke = 0.25,     # 黑色边线粗细      ← 修改此处（0 = 无边线）
    color  = "black", # 边线颜色固定为黑色 ← 修改此处换其他边线颜色
    alpha  = 0.9
  ) +
  scale_fill_manual(values = GROUP_COLORS[setdiff(GROUP_PREFIXES, CONTROL_GROUP)]) +
  labs(x = NULL, y = NULL, title = "Log2 Fold Change") +
  theme_bw() +
  theme(
    axis.text.y      = element_text(face = "italic", size = 10),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 11, face = "bold")
  )


# ── P2：右图 — 相对丰度水平柱状图 ──────────────────────────────────
# 要求：
#   - 删除背景底纹（panel.grid 全部去掉）
#   - 只保留外边框（panel.border 保留，panel.background 透明）
#   - 误差棒粗细适中
#   - 组顺序：PTW（红）→ 1 dpi → 3 dpi → 5 dpi，从上到下排列

p2 <- ggplot(
    plot_data,
    aes(x = Taxon, y = Mean, fill = Group, group = Group)
  ) +
  # 柱子：reverse=TRUE 使得 coord_flip 后从上到下为 PTW/1 dpi/3 dpi/5 dpi
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = DODGE_WIDTH, reverse = TRUE),
    width    = BAR_WIDTH
  ) +
  # 误差棒：与柱子 dodge 参数完全一致
  geom_errorbar(
    aes(ymin = ifelse(Mean - SD < 0, 0, Mean - SD),
        ymax = Mean + SD),
    position  = position_dodge(width = DODGE_WIDTH, reverse = TRUE),
    width     = ERRORBAR_CAPWIDTH,    # ← 帽子宽度，在参数区调整
    linewidth = ERRORBAR_LINEWIDTH,   # ← 线条粗细，在参数区调整
    color     = "black"
  ) +
  coord_flip() +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = NULL, title = "RA (%)\nmean ± SD") +
  theme_bw() +
  theme(
    # 隐藏 Y 轴文字和刻度（左图已有属名）
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    # 删除所有背景网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 只保留外边框（panel.border 默认已有，不用额外设置）
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, size = 11, face = "bold")
  )


# ── 提取公共图例（从 P2）────────────────────────────────────────────
shared_legend <- get_legend(
  p2 + theme(
    legend.position    = "bottom",
    legend.justification = "center",
    legend.direction   = "horizontal",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.5, "cm")
  )
)


# ── 组合图片 ────────────────────────────────────────────────────────
# rel_widths：左图需要更多宽度（含属名），右图较窄
p_combined <- plot_grid(
  p1, p2,
  ncol       = 2,
  rel_widths = c(1.5, 0.9),
  align      = "h",
  axis       = "bt"
)

final_plot <- plot_grid(
  p_combined, shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.08)
)

final_plot
# ── 输出 ────────────────────────────────────────────────────────────
print(final_plot)

if (!is.null(OUTPUT_PDF)) {
  ggsave(OUTPUT_PDF, final_plot,
         width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT)
  message("图片已保存: ", OUTPUT_PDF)
}

if (!is.null(OUTPUT_RA_TABLE)) {
  write.table(
    ra_percent_table,
    file = OUTPUT_RA_TABLE,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("RA百分比统计表已保存: ", OUTPUT_RA_TABLE)
}

# 导出宽格式表（每个 Group 为列，保留 mean 和 sd 的数值列）
ra_percent_table_wide <- ra_percent_table %>%
  select(Taxon, Group, RA_percent_mean, RA_percent_sd) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(RA_percent_mean, RA_percent_sd),
    names_sep = "_"
  ) %>%
  arrange(factor(Taxon, levels = taxon_order))

if (!is.null(OUTPUT_RA_TABLE_WIDE)) {
  write.table(
    ra_percent_table_wide,
    file = OUTPUT_RA_TABLE_WIDE,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("宽格式RA百分比统计表已保存: ", OUTPUT_RA_TABLE_WIDE)
}

# 导出格式化（mean ± SD）表，便于直接查看/报告
ra_percent_table_formatted_wide <- ra_percent_table %>%
  mutate(fmt = sprintf("%.4f%% ± %.4f%%", RA_percent_mean, RA_percent_sd)) %>%
  select(Taxon, Group, fmt) %>%
  pivot_wider(names_from = Group, values_from = fmt) %>%
  arrange(factor(Taxon, levels = taxon_order))

if (!is.null(OUTPUT_RA_TABLE_FORMATTED)) {
  write.table(
    ra_percent_table_formatted_wide,
    file = OUTPUT_RA_TABLE_FORMATTED,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("格式化RA百分比表已保存: ", OUTPUT_RA_TABLE_FORMATTED)
}
