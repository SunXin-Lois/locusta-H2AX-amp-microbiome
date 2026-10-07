## This script runs on R version 4.2 or higher

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggstance)
library(grid)

# 请从本脚本所在的“绘图数据”目录运行。

# 1. 读取更新后的数据文件（不修改原始文件）
#tree_file <- "ASV_LmH_g_Unclassified_f_Enterobacteriaceae_filtered_PTW3.nwk"
tree_file <- "ASV_LmH_g_Unclassified_f_Enterobacteriaceae_filtered_CK3_2.nwk"
table_file <- "table_filtered_w_LmH_Unclassified_f_Enterobacteriaceae_updated_PTW_filtered3.tsv"
output_file <- "NJ_tree_Unclassified_f_Enterobacteriaceae_RA_PTW_1dpi_3dpi_5dpi_20260901.pdf"

tree <- read.tree(tree_file)
df_raw <- read.table(
  table_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

# 参考图将 0–1 的节点支持度显示为整数（例如 0.8410 显示为 841）
tree$node.label <- as.character(round(as.numeric(tree$node.label) * 1000))

# 2. 计算全局总丰度（所有 ASV 在所有样本中的总和）
total_abundance_all <- sum(df_raw)

# 3. 计算各组平均值
group_cols <- list(
  PTW = c("PTW-1", "PTW-2", "PTW-3", "PTW-4", "PTW-5"),
  `1 dpi` = c("1 dpi-1", "1 dpi-2", "1 dpi-3", "1 dpi-4", "1 dpi-5"),
  `3 dpi` = c("3 dpi-1", "3 dpi-2", "3 dpi-3", "3 dpi-4", "3 dpi-5"),
  `5 dpi` = c("5 dpi-1", "5 dpi-2", "5 dpi-3", "5 dpi-4", "5 dpi-5")
)

mean_df <- data.frame(ID = rownames(df_raw))
for (g_name in names(group_cols)) {
  mean_df[[g_name]] <- rowMeans(df_raw[, group_cols[[g_name]]])
}

# 4. 转换为长格式并计算占全体总丰度的百分比
plot_data <- mean_df %>%
  pivot_longer(cols = -ID, names_to = "Group", values_to = "MeanCount") %>%
  mutate(Global_RA = (MeanCount / total_abundance_all) * 100)

# 设置堆叠顺序：从左到右依次为 PTW、1 dpi、3 dpi、5 dpi
plot_data$Group <- factor(
  plot_data$Group,
  levels = rev(c("PTW", "1 dpi", "3 dpi", "5 dpi"))
)

# 5. 绘图
tree_width <- max(fortify(tree)$x)
highlight_tips <- c("ASV_4260", "ASV_7121")
panel_name <- "Re. Abundance (%)"

p <- ggtree(tree, color = "black") +
  geom_text2(
    aes(subset = !isTip, label = label),
    family = "Arial",
    size = 1.1,
    hjust = 1.15,
    vjust = -0.45,
    color = "#222222"
  ) +
  geom_tiplab(
    aes(color = label %in% highlight_tips),
    align = TRUE,
    linetype = "dotted",
    linesize = 0.3,
    family = "Arial",
    size = 2.55,
    offset = tree_width * 0.1
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "#ED1C24"),
    guide = "none"
  )

p_final <- facet_plot(
  p,
  panel = panel_name,
  data = plot_data,
  geom = geom_barh,
  mapping = aes(x = Global_RA, fill = Group),
  stat = "identity",
  width = 0.7
) +
  scale_fill_manual(
    values = c(
      "PTW" = "#999999",
      "1 dpi" = "#E69F00",
      "3 dpi" = "#56B4E9",
      "5 dpi" = "#D55E00"
    ),
    breaks = c("PTW", "1 dpi", "3 dpi", "5 dpi"),
    name = "Group"
  ) +
  scale_x_continuous(
    position = "top",
    breaks = c(0, 10, 20),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  xlim_tree(tree_width * 4) +
  theme_tree2() +
  theme(
    text = element_text(family = "Arial", color = "#050504"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(8, "pt"),
    legend.spacing.y = unit(1, "pt"),
    legend.box.spacing = unit(4, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.5, face = "bold"),
    axis.text.x.top = element_text(size = 7.5, color = "#222222"),
    axis.ticks.x.top = element_line(linewidth = 0.3, color = "#222222"),
    axis.line.x.top = element_line(linewidth = 0.3, color = "#222222"),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )

# 参考图中树面板和丰度面板的宽度比约为 1.9:1
p_final <- facet_widths(
  p_final,
  widths = c(Tree = 1.9, "Re. Abundance (%)" = 1)
)

# 6. 输出 PDF。参考图页面为 171.311 × 147.649 pt
# （即约 2.379 × 2.051 英寸）。
print(p_final)

figure_width <- 171.311 / 72
figure_height <- 147.649 / 72

if (capabilities("cairo")) {
  cairo_pdf(
    filename = output_file,
   # width = figure_width,
   # height = figure_height,
    family = "Arial",
    bg = "white"
  )
} else {
  pdf(
    file = output_file,
    #width = figure_width,
    #height = figure_height,
    family = "sans",
    bg = "white"
  )
}
print(p_final)
dev.off()
