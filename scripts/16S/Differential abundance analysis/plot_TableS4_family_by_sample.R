## Table S4 科水平相对丰度 — 每个样本单独展示的堆叠条形图
## 数据：TableS4_Relative_abundance_Fam.xls

.libPaths(c(.libPaths(), Sys.getenv("R_LIBS_USER")))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ── 参数 ────────────────────────────────────────────────────────────
INPUT_FILE <- "TableS4_Relative_abundance_Fam.xls"

# 堆叠顺序：从下到上
TAXON_ORDER <- c(
  "Enterobacteriaceae",
  "Streptococcaceae",
  "Leuconostocaceae",
  "Enterococcaceae",
  "Yersiniaceae",
  "Unclassified",
  "Others"
)

TAXON_COLORS <- c(
  "Enterobacteriaceae" = "#E8735A",
  "Streptococcaceae"   = "#39B5B2",
  "Leuconostocaceae"   = "#1A7A78",
  "Enterococcaceae"    = "#7B6FAB",
  "Yersiniaceae"       = "#5A8FC4",
  "Unclassified"       = "#8E8E8E",
  "Others"             = "#D4CEBC"
)

BAR_WIDTH <- 0.85
FIG_HEIGHT <- 4.2
DPI_OUT <- 300


# ── 读取并整理 ──────────────────────────────────────────────────────
raw_df <- read_excel(INPUT_FILE, skip = 1)
names(raw_df)[1] <- "Sample"

plot_long <- raw_df %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Group = case_when(
      grepl("^PTW-", Sample) ~ "PTW",
      grepl("^1 dpi-", Sample) ~ "1 dpi",
      grepl("^3 dpi-", Sample) ~ "3 dpi",
      grepl("^5 dpi-", Sample) ~ "5 dpi",
      grepl("^RNAi-", Sample) ~ "RNAi",
      grepl("^Rescue-", Sample) ~ "Rescue",
      TRUE ~ NA_character_
    ),
    Replicate = factor(
      sub(".*-", "", Sample),
      levels = as.character(1:5)
    )
  ) %>%
  filter(!is.na(Group))

plot_long$Taxon <- factor(plot_long$Taxon, levels = TAXON_ORDER)


# ── 图 1：PTW / 1 dpi / 3 dpi / 5 dpi ───────────────────────────────
plot_dpi <- plot_long %>%
  filter(Group %in% c("PTW", "1 dpi", "3 dpi", "5 dpi")) %>%
  mutate(Group = factor(Group, levels = c("PTW", "1 dpi", "3 dpi", "5 dpi")))

p_dpi <- ggplot(plot_dpi, aes(x = Replicate, y = Abundance, fill = Taxon)) +
  geom_col(
    width = BAR_WIDTH,
    color = "white",
    linewidth = 0.2,
    # ggplot2 默认堆叠方向与参考图相反，需 reverse = TRUE
    # 这样 Enterobacteriaceae 才会出现在柱子底部
    position = position_stack(reverse = TRUE)
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_manual(values = TAXON_COLORS, breaks = TAXON_ORDER) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Replicate", y = "Relative abundance", fill = "Family") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  )

print(p_dpi)

ggsave(
  "Family_abundance_S4_by_sample_PTW_1_3_5dpi.pdf",
  plot = p_dpi,
  width = 10.5,
  height = FIG_HEIGHT,
  bg = "white"
)
ggsave(
  "Family_abundance_S4_by_sample_PTW_1_3_5dpi_fixed.png",
  plot = p_dpi,
  width = 10.5,
  height = FIG_HEIGHT,
  dpi = DPI_OUT,
  bg = "white"
)

write.table(
  plot_dpi,
  "Family_abundance_S4_by_sample_PTW_1_3_5dpi_long.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mean_sd_dpi <- plot_dpi %>%
  group_by(Group, Taxon) %>%
  summarise(
    N = n_distinct(Sample),
    RA_percent_mean = mean(Abundance, na.rm = TRUE) * 100,
    RA_percent_sd = sd(Abundance, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(factor(Group, levels = c("PTW", "1 dpi", "3 dpi", "5 dpi")),
          factor(Taxon, levels = TAXON_ORDER))

write.table(
  mean_sd_dpi,
  "Family_abundance_S4_by_sample_PTW_1_3_5dpi_mean_SD_by_group.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("图 1 已保存: Family_abundance_S4_by_sample_PTW_1_3_5dpi_fixed.png / .pdf")


# ── 图 2：3 dpi / RNAi / Rescue ─────────────────────────────────────
plot_fig2 <- plot_long %>%
  filter(Group %in% c("3 dpi", "RNAi", "Rescue")) %>%
  mutate(Group = factor(Group, levels = c("3 dpi", "RNAi", "Rescue")))

p_fig2 <- ggplot(plot_fig2, aes(x = Replicate, y = Abundance, fill = Taxon)) +
  geom_col(
    width = BAR_WIDTH,
    color = "white",
    linewidth = 0.2,
    position = position_stack(reverse = TRUE)
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_manual(values = TAXON_COLORS, breaks = TAXON_ORDER) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Replicate", y = "Relative abundance", fill = "Family") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  )

print(p_fig2)

ggsave(
  "Family_abundance_S4_by_sample_3dpi_RNAi_Rescue.pdf",
  plot = p_fig2,
  width = 8.5,
  height = FIG_HEIGHT,
  bg = "white"
)
ggsave(
  "Family_abundance_S4_by_sample_3dpi_RNAi_Rescue.png",
  plot = p_fig2,
  width = 8.5,
  height = FIG_HEIGHT,
  dpi = DPI_OUT,
  bg = "white"
)

write.table(
  plot_fig2,
  "Family_abundance_S4_by_sample_3dpi_RNAi_Rescue_long.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mean_sd_fig2 <- plot_fig2 %>%
  group_by(Group, Taxon) %>%
  summarise(
    N = n_distinct(Sample),
    RA_percent_mean = mean(Abundance, na.rm = TRUE) * 100,
    RA_percent_sd = sd(Abundance, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(factor(Group, levels = c("3 dpi", "RNAi", "Rescue")),
          factor(Taxon, levels = TAXON_ORDER))

write.table(
  mean_sd_fig2,
  "Family_abundance_S4_by_sample_3dpi_RNAi_Rescue_mean_SD_by_group.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 导出宽格式原始表（供其他脚本复用）
write.table(
  raw_df,
  "Table_S4_family_relative_abundance.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("图 2 已保存: Family_abundance_S4_by_sample_3dpi_RNAi_Rescue.png / .pdf")
message("配套表已保存: *_long.tsv, *_mean_SD_by_group.tsv, Table_S4_family_relative_abundance.tsv")
