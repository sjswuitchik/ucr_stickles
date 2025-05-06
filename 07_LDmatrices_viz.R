library(reshape2)
library(viridis)
library(tidyverse)

mce_ld_mat <- read.table("mce_ld_matrix.ld", header = F)
rownames(mce_ld_mat) <- mce_sig_ann$ID
colnames(mce_ld_mat) <- mce_sig_ann$ID

mce_ld_long <- melt(as.matrix(mce_ld_mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

ggplot(mce_ld_long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²)", x = "SNP", y = "SNP", fill = expression(r^2))

md_ld_mat <- read.table("md_ld_matrix.ld", header = F)
rownames(md_ld_mat) <- md_sig_ann$ID
colnames(md_ld_mat) <- md_sig_ann$ID

md_ld_long <- melt(as.matrix(md_ld_mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

ggplot(md_ld_long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²)", x = "SNP", y = "SNP", fill = expression(r^2))
