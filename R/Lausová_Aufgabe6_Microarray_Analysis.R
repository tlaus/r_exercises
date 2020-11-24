## Aufgabe 6 - Microarray analysis
## Tereza Lausova
## 15.9.2020

# Load and setup libraries
library(here)
library(tidyverse)
# install.packages("O:/Lets_do_research_on_corona/hgu133plus2hsenstcdf_24.0.0.tar.gz",repos = NULL)
# install.packages("O:/Lets_do_research_on_corona/hgu133plus2hsenstprobe_24.0.0.tar.gz",repos = NULL)
library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)
# BiocManager::install("vsn")
library(vsn)
# browseVignettes("vsn")
# BiocManager::install("affy")
library(affy)
library(AnnotationDbi)
library(reshape2)
library(biomaRt)
library(qpdf)

# load data
# bc_cel <- ReadAffy(celfile.path = here("rawdata"))
# bc_cel@cdfName <- "HGU133Plus2_Hs_ENST"

# save the data as RDS
# saveRDS(bc_cel, here("rds/bc_cel.RDS"))

# QC
bc_cel <- readRDS(here("rds/bc_cel.RDS"))

## export all images to a pdf for inspection
image(bc_cel, col = rainbow(100, start=0, end=0.75)[100:1])

# GSM687020 - green flake on the scan
singlechip <- bc_cel[5]
pdf(here("plots/Lausova_Aufgabe6_QC_singlechip.pdf"))
image(singlechip, col = rainbow(100, start=0, end=0.75)[100:1])
dev.off()

# pheno data
samples <- pData(bc_cel) %>%
  rownames() %>%
  str_sub(1,9)

# normalization
bc_norm <- vsnrma(bc_cel)
saveRDS(bc_norm, here("rds/bc_norm.RDS"))

## verify fit
meanSdPlot(bc_norm)

## before, after normalization comparison
str(exprs(bc_cel))
class(exprs(bc_cel))

### before
tidy_exprs <- exprs(bc_cel) %>%
  melt() %>%
  as_tibble() %>%
  mutate("sample" = Var2)

pdf(here("plots/Lausova_Aufgabe6_rawdata.pdf"))
ggplot(tidy_exprs, aes(x = log(value), fill = sample, y = sample)) +
  theme_minimal() +
  geom_boxplot() +
  guides(fill = "none") +
  labs(title = "Breast cancer microarray data before normalization")
  # facet_wrap(~sample, nrow = 10) # for density plots
dev.off()

### after
tidy_norm <- exprs(bc_norm) %>%
  melt() %>%
  as_tibble() %>%
  mutate("sample" = Var2)

pdf(here("plots/Lausova_Aufgabe6_normalized.pdf"))
ggplot(tidy_norm, aes(x = log(value), fill = sample, y = sample)) +
  theme_minimal() +
  geom_boxplot() +
  guides(fill = "none")+
  labs(title = "Breast cancer microarray data after normalization")
dev.off()

# some are sligthly bimodal even after normalization

## RNA degradation plot
bc_rnadeg <- AffyRNAdeg(bc_cel)
pdf(here("plots/Lausova_Aufgabe6_RNAdeg.pdf"))
plotAffyRNAdeg(bc_rnadeg, cols = rainbow(10))
plotAffyRNAdeg(bc_rnadeg, cols = rainbow(10), transform = "shift.only")
dev.off()

## scatterplot
pdf(here("plots/Lausova_Aufgabe6_QC_scatterplot.pdf"))
map(c(1:length(colnames(bc_norm))-1), function(sample) {
  plot(exprs(bc_norm)[,c(sample, sample+1)], pch=".")
  abline(0,1,col="red")
})
dev.off()
length_pdf <- pdf_length(input = here("plots/Lausova_Aufgabe6_QC_scatterplot.pdf"))
pdf_subset(input = here("plots/Lausova_Aufgabe6_QC_scatterplot.pdf"), pages = 2:length_pdf, output = NULL)

## density plots

pdf(here("plots/Lausova_Aufgabe6_histogram.pdf"), onefile = T, width = 6, height = 4)
## Not normalized
ggplot(tidy_exprs) +
  theme_minimal() +
  geom_density(aes(x = value), fill = "cornflowerblue") +
  labs(title = "Density plot of non-normalised expression data")

## Normalized
ggplot(tidy_norm) +
  theme_minimal() +
  geom_density(aes(x = value), fill = "cornflowerblue")+
  labs(title = "Density plot of normalised expression data")

dev.off()

# analysis
## add annotation
rownames(bc_norm) %>% View() # first 62 are affymetrix

tr_names <- rownames(bc_norm)[63:length(rownames(bc_norm))] %>%
  str_replace("_at", "")
all_names <- rownames(bc_norm) %>%
  str_replace("_at", "")

ensembl <- useMart("ensembl")
hs <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
entrez_data <- getBM(mart = hs, 
                     attributes = c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "start_position", "entrezgene_id", "external_gene_name", "external_transcript_name", "family", "ensembl_transcript_id_version"), 
                     filters = c("ensembl_transcript_id_version"), values = tr_names)

length(tr_names)
length(all_names)
length(entrez_data$ensembl_transcript_id)
bc_filtered <- exprs(bc_norm)[which(all_names %in% entrez_data$ensembl_transcript_id_version), ]
all_names <- all_names[which(all_names %in% entrez_data$ensembl_transcript_id_version)]
rownames(bc_filtered) <- entrez_data$ensembl_transcript_id_version[which(all_names %in% entrez_data$ensembl_transcript_id_version)]

# matched <- which(entrez_data$ensembl_transcript_id_version %in% rownames(bc_filtered)) # This is not correct (probably)
# bc_norm_annotated <- cbind(entrez_data[matched,], bc_filtered)
# write_csv(entrez_data, here("tables/Lausova_Aufgabe6_Microarray_ILgenes.csv"))

## gene expression
entrez_data_filtered <- entrez_data[which(all_names %in% entrez_data$ensembl_transcript_id_version),]
rownames(bc_filtered) <- entrez_data_filtered$external_transcript_name
interleukins_boolean <- grepl("^IL\\d+", rownames(bc_filtered))
any(interleukins_boolean)
bc_IL <- bc_filtered[interleukins_boolean, ]

bc_IL_annotated <- cbind(bc_IL, entrez_data_filtered[interleukins_boolean,])
write_csv(bc_IL_annotated, "tables/Lausova_Aufgabe6_Microarray_ILgenes.csv")

bc_all_annotated <- cbind(bc_filtered, entrez_data_filtered)
write_csv(bc_all_annotated, "tables/Lausova_Aufgabe6_Microarray_Genes.csv")


bc_IL_tidy <- bc_IL %>%
  melt() %>%
  as_tibble() %>%
  mutate("transcript" = as.character(Var1), "sample" = Var2, "intensity" = value) %>%
  # arrange(desc(transcript)) %>%
  mutate(transcript_groups = str_extract(transcript, pattern = "[:alnum:]{2,}-") %>% str_sub(end = -2)) %>%
  arrange(transcript_groups) %>%
  group_by(transcript_groups) %>%
  arrange(transcript) %>%
  ungroup() %>%
  # group_by(transcript) %>%
  # mutate(median_val = median(value)) %>%
  # arrange(desc(median_val)) %>%
  # ungroup() %>%
  mutate_if(is.character, as.factor)


pdf(here("plots/Lausova_Aufgabe6_Gene_expression_IL_breast_cancer.pdf"), height = 36, width = 6)
ggplot(bc_IL_tidy, aes(x = intensity, y = transcript, fill = transcript_groups)) +
  geom_boxplot() +
  guides(fill = "none")
dev.off()

# heatmap
library(pheatmap)

pheatmap(bc_IL)
