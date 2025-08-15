#### Create LD matrices that are plotable

mce.ld.mat <- read.table("mce_ld_matrix.ld", header = F)
rownames(mce.ld.mat) <- mce.sig.ann$variant.id
colnames(mce.ld.mat) <- mce.sig.ann$variant.id

mce.ld.long <- melt(as.matrix(mce.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

md.ld.mat <- read.table("md_ld_matrix.ld", header = F)
rownames(md.ld.mat) <- md.sig.ann$variant.id
colnames(md.ld.mat) <- md.sig.ann$variant.id

md.ld.long <- melt(as.matrix(md.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

time_hdvmg.ld.mat <- read.table("time_hdvmg_ld_matrix.ld", header = F)
rownames(time_hdvmg.ld.mat) <- time.hdvmg.sig.ann$variant.id
colnames(time_hdvmg.ld.mat) <- time.hdvmg.sig.ann$variant.id

time_hdvmg.ld.long <- melt(as.matrix(time_hdvmg.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

time_mdvmg.ld.mat <- read.table("time_mdvmg_ld_matrix.ld", header = F)
rownames(time_mdvmg.ld.mat) <- time.mdvmg.sig.ann$variant.id
colnames(time_mdvmg.ld.mat) <- time.mdvmg.sig.ann$variant.id

time_mdvmg.ld.long <- melt(as.matrix(time_mdvmg.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

ttpg.ld.mat <- read.table("ttpg_ld_matrix.ld", header = F)
rownames(ttpg.ld.mat) <- ttpg.sig.ann$variant.id
colnames(ttpg.ld.mat) <- ttpg.sig.ann$variant.id

ttpg.ld.long <- melt(as.matrix(ttpg.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

## Heat maps
ggplot(mce.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) MaxCranElev", x = "SNP", y = "SNP", fill = expression(r^2))

ggplot(md.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) MaxDecel", x = "SNP", y = "SNP", fill = expression(r^2))

ggplot(time_hdvmg.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) Time HDvMG", x = "SNP", y = "SNP", fill = expression(r^2))

ggplot(time_mdvmg.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) Time MDvMG", x = "SNP", y = "SNP", fill = expression(r^2))

ggplot(ttpg.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) TTPG", x = "SNP", y = "SNP", fill = expression(r^2))


##########################
## Morphological Traits ##
##########################

ray.ld.mat <- read.table("gwas_out/ray_ld_matrix.ld", header = F)
rownames(ray.ld.mat) <- ray.sig.ann$variant.id
colnames(ray.ld.mat) <- ray.sig.ann$variant.id

ray.ld.long <- melt(as.matrix(ray.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

ggplot(ray.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) Ray Count", x = "SNP", y = "SNP", fill = expression(r^2))
