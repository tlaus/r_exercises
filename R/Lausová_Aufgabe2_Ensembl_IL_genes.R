## Aufgabe 2 - Ensembl IL genes
## Tereza Lausová
## 14.9.2020

# Load libraries
library(tidyverse)
library(here)
if (!require(biomaRt)) {
  BiocManager::install("biomaRt")
}

# Connect to the Ensembl
ensembl <- useMart("ensembl")
hs <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

# listFilters(hs)
# listAttributes(hs)

# IL gene family obtained from https://www.genenames.org/data/genegroup/#!/group/601
interleukins <- read_csv(here("data/group-601.csv"), skip = 1)
# interleukins$`Approved symbol`

# Obtain IL info from Ensembl
hs_data <- getBM(mart = hs, 
                 attributes = c("ensembl_gene_id", "chromosome_name", "start_position","external_gene_name"), 
                 filters = c("external_gene_name"), 
                 values = interleukins$`Approved symbol`)

# Assign chromosome names as factors
hs_data$chromosome_name <- factor(hs_data$chromosome_name, levels = c(1:22, "X", "Y"))
hs_data$external_gene_name_levels <- as.numeric(as.factor(hs_data$external_gene_name))

# Plot 
pdf(file = here("plots/Lausová_Aufgabe2_Ensembl_IL_genes.pdf"))
ggplot(data = hs_data, aes(x = chromosome_name)) + 
  theme_minimal() +
  geom_bar(aes(fill = chromosome_name), stat = "count", alpha = 0.5,
           position = ) +
  geom_point(aes(y = external_gene_name_levels/4)) +
  geom_segment(aes(y = external_gene_name_levels/4, 
                   yend = external_gene_name_levels/4,
                   x = chromosome_name,
                   xend = 25, 
                   group = external_gene_name), linetype = "dashed") +
  scale_x_discrete(drop = F) +
  scale_y_continuous(breaks = 1:(max(summary.factor(hs_data$chromosome_name))*4),
                     sec.axis = sec_axis(~.*4,
                                         breaks = 1:length(hs_data$external_gene_name_levels),
                                         labels = hs_data$external_gene_name)) +
  guides(fill = "none") +
  labs(title = "Interleukin family genes location on chromosomes")
dev.off()

# Table export
write_csv(hs_data, path = here("tables/Lausová_Aufgabe2_Ensembl_IL_genes.csv"))