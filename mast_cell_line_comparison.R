# Dan Dwyer paper on subsets of mast cells. how similar is our mast cell line cell to the subtypes.
# https://pubmed.ncbi.nlm.nih.gov/39744949/

library(tidyverse)
# library(BIGpicture)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(edgeR)
library(limma)
library(ggpubr)
# library(RNAetc)

#### Gene key ####
ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                               dataset="hsapiens_gene_ensembl")
key <- biomaRt::getBM(attributes=c("ensembl_gene_id", "hgnc_symbol",
                                   "gene_biotype"),
                      mart=ensembl)

#### Load and format P546-1 data ####
load("dat_clean/voom_objects/mc_voom.RData")

#Add gene metadata
# load("dat_clean/ensembl_gene_dictionary.RData")
mc_voom$genes <- key %>% 
  filter(ensembl_gene_id %in% rownames(mc_voom$E)) %>% 
  mutate(ensembl_gene_id = factor(
    ensembl_gene_id, 
    levels = rownames(mc_voom$E))) %>%
  arrange(ensembl_gene_id)

#Change to hgnc symbol
# mc_voom_hgnc <- mc_voom
# mc_voom_hgnc$E <- as.data.frame(mc_voom_hgnc$E) %>% 
#   rownames_to_column("ensembl_gene_id") %>% 
#   left_join(ens_gene_dict %>% select(ensembl_gene_id,geneName)) %>% 
#   select(-ensembl_gene_id) 

#### Load and format Dwyer data ####
#pseudo bulk counts
mast <- read_excel("dat_other/Dwyer JCI MC Subset Pseudobulk.xlsx",
                   sheet=1)
mast <- mast[,c(2,9:ncol(mast))]

#create metadata
mast_meta <- data.frame(libID=colnames(mast)[-1]) %>% 
  separate_wider_delim(libID, names=c("Stimulation","Donor.ID"), 
                       delim="_Donor", cols_remove = FALSE) %>% 
  mutate(Treatment="Alone",Timepoint="Baseline", study="Dwyer")

#Rename genes to ensembl
mast_ens <- mast %>% 
  inner_join(key, by=c("Gene"="hgnc_symbol")) %>% 
  select(-Gene, -gene_biotype) %>% 
  select(ensembl_gene_id, everything()) %>% 
  column_to_rownames("ensembl_gene_id")

mast_genes <- key %>% 
  filter(ensembl_gene_id %in% rownames(mast_ens)) %>% 
  mutate(ensembl_gene_id = factor(
    ensembl_gene_id, 
    levels = rownames(mast_ens))) %>%
  arrange(ensembl_gene_id)

#Limma normalize
mast_dat <- edgeR::DGEList(
  counts = mast_ens,
  samples = mast_meta, 
  genes = mast_genes)
mast_dat_norm <- calcNormFactors(mast_dat)
mast_dat_voom <- voomWithQualityWeights(
  mast_dat_norm,
  design = model.matrix(~Stimulation, mast_dat_norm$samples))
rm(mast, mast_ens, mast_meta, mast_genes, mast_dat, mast_dat_norm)

#### Filter key ####
key_filter <- key %>% 
  filter(hgnc_symbol %in% mc_voom$genes$hgnc_symbol & 
           hgnc_symbol %in% mast_dat_voom$genes$hgnc_symbol)

#### Marker genes ####
markers <- data.frame()
sheet_names <- excel_sheets("dat_other/Dwyer JCI MC Subset Pseudobulk.xlsx")
for(i in sheet_names){
  markers <- read_excel("dat_other/Dwyer JCI MC Subset Pseudobulk.xlsx",
                        sheet=i) %>% 
    select(Gene:padj) %>% 
    mutate(group = i) %>% 
    bind_rows(markers)
}

markers_format <- markers %>% 
  mutate(group = gsub("DESeq2-pseudobulk-","",group),
         group = gsub(" vs all","",group)) %>% 
  rename(hgnc_symbol=Gene) 

markers_top50 <- markers_format  %>% 
  filter(log2FoldChange>0) %>% 
  slice_head(n=50, by=group) %>% 
  pull(hgnc_symbol) %>% unique()
markers_top100 <- markers_format %>% 
  filter(log2FoldChange>0) %>% 
  slice_head(n=100, by=group) %>% 
  pull(hgnc_symbol) %>% unique()
markers_top200 <- markers_format %>% 
  filter(log2FoldChange>0) %>% 
  slice_head(n=200, by=group) %>% 
  pull(hgnc_symbol) %>% unique()

#convert to ens
markers_top50_ens <- key_filter %>% 
  filter(hgnc_symbol %in% markers_top50) %>% 
  filter(ensembl_gene_id %in% mc_voom$genes$ensembl_gene_id & 
           ensembl_gene_id %in% mast_dat_voom$genes$ensembl_gene_id) %>% 
  pull(ensembl_gene_id)
markers_top100_ens <- key_filter %>% 
  filter(hgnc_symbol %in% markers_top100) %>% 
  filter(ensembl_gene_id %in% mc_voom$genes$ensembl_gene_id & 
           ensembl_gene_id %in% mast_dat_voom$genes$ensembl_gene_id) %>% 
  pull(ensembl_gene_id)

markers_top200_ens <- key_filter %>% 
  filter(hgnc_symbol %in% markers_top200) %>% 
  filter(ensembl_gene_id %in% mc_voom$genes$ensembl_gene_id & 
           ensembl_gene_id %in% mast_dat_voom$genes$ensembl_gene_id) %>% 
  pull(ensembl_gene_id)

#### Average pseudobulk ####
mast_dat_voom_ave <- as.data.frame(mast_dat_voom$E) %>%
  rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id, names_to = "libID") %>%
  separate_wider_delim(libID, names=c("Stimulation","Donor.ID"),
                       delim="_Donor", cols_remove = FALSE) %>%
  group_by(Stimulation,ensembl_gene_id) %>%
  summarise(value=mean(value, na.rm=TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Stimulation) %>%
  column_to_rownames("ensembl_gene_id")

#### Correlation ####
dat_markers_mc50 <- mc_voom$E[markers_top50_ens,]
dat_markers_mast50 <- mast_dat_voom_ave[markers_top50_ens,]
dat_markers_mc100 <- mc_voom$E[markers_top100_ens,]
dat_markers_mast100 <- mast_dat_voom_ave[markers_top100_ens,]
dat_markers_mc200 <- mc_voom$E[markers_top200_ens,]
dat_markers_mast200 <- mast_dat_voom_ave[markers_top200_ens,]

sim_matrix50 <- cor(dat_markers_mast50, dat_markers_mc50,
                    method = "pearson")
sim_matrix100 <- cor(dat_markers_mast100, dat_markers_mc100,
                    method = "pearson")
sim_matrix200 <- cor(dat_markers_mast200, dat_markers_mc200,
                     method = "pearson")

#### Assign type ####
# assign each bulk sample to closest pseudobulk cluster
assigned_cluster50 <- apply(sim_matrix50, 2, 
                          function(x) rownames(sim_matrix50)[which.max(x)])
assigned_cluster100 <- apply(sim_matrix100, 2, 
                            function(x) rownames(sim_matrix100)[which.max(x)])
assigned_cluster200 <- apply(sim_matrix200, 2, 
                             function(x) rownames(sim_matrix200)[which.max(x)])

sim_result <- mc_voom$targets %>% 
  left_join(data.frame(assigned_cluster50) %>%
              rownames_to_column("libID")) %>% 
  left_join(data.frame(assigned_cluster100) %>%
             rownames_to_column("libID")) %>% 
  left_join(data.frame(assigned_cluster200) %>%
              rownames_to_column("libID")) %>% 
  select(libID, Stimulation, Timepoint, Treatment, 
         assigned_cluster50, assigned_cluster100,
         assigned_cluster200) %>% 
  separate_wider_delim(assigned_cluster50, names="cluster50",
                       delim="_Donor", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  separate_wider_delim(assigned_cluster100, names="cluster100",
                       delim="_Donor", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  separate_wider_delim(assigned_cluster200, names="cluster200",
                       delim="_Donor", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  rowwise() %>% 
  mutate(concordance = ifelse(cluster50==cluster100 & 
                                cluster100==cluster200, TRUE, FALSE))
View(sim_result)

#### Summary table ####
as.data.frame(sim_matrix100) %>% 
  rownames_to_column() %>% 
  separate_wider_delim(rowname, 
                       names=c("cluster"), 
                       too_many = "drop",
                       delim="_Donor", cols_remove = FALSE) %>% 
  pivot_longer(-c(cluster,rowname), names_to = "libID") %>% 
  group_by(libID, cluster) %>% 
  summarise(meanCor = signif(mean(value), digits=2),
            sdCor = signif(sd(value), digits=2),
            cluster_summ = paste(meanCor,sdCor,sep="+")) %>% 
  ungroup() %>% 
  left_join(sim_result) %>% 
  select(libID, Stimulation, Timepoint,
         Treatment, cluster, cluster_summ, cluster100) %>% 
  pivot_wider(names_from = cluster, values_from = cluster_summ) %>% view

#### Heatmap ####
hm_dat <- as.data.frame(sim_matrix100) %>% 
  rownames_to_column() %>% 
  separate_wider_delim(rowname, 
                       names=c("cluster"), 
                       too_many = "drop",
                       delim="_Donor", cols_remove = FALSE) %>% 
  pivot_longer(-c(cluster,rowname), names_to = "libID") %>% 
  group_by(libID, cluster) %>% 
  summarise(meanCor = mean(value)) %>% 
  ungroup() %>% 
  left_join(sim_result) %>% 
  select(libID, Stimulation, Timepoint,
         Treatment, cluster, meanCor, cluster100) %>% 
  pivot_wider(names_from = cluster, values_from = meanCor) %>% 
  column_to_rownames("libID") %>% 
  arrange(Stimulation, Treatment, Timepoint) %>% 
  mutate(Stimulation = recode(Stimulation, "HRV Infection"="HRV"),
         Treatment = recode(Treatment, "AECs"="AEC coculture"),
         Timepoint = recode(Timepoint, "Baseline"="0"),
         Timepoint = gsub(" hrs","",Timepoint),
         Timepoint = factor(Timepoint, levels=c("0","4","8",
                                                "28","52")))

hm_mat <- as.matrix(hm_dat[,5:10])
# Sample annotation
# sample_annot <- HeatmapAnnotation(
#   Infection = hm_dat$Stimulation,
#   Coculture = hm_dat$Treatment,
#   Hour = hm_dat$Hour,
#   col = list(Infection = c("None" = "skyblue", "HRV" = "salmon"),
#              Coculture = c("Alone"="red", "AECs"="black"),
#              Hour = c("0"="grey90",  "4"="grey70",  "8"="grey50",
#                       "28"="grey30", "52"="grey10"))
# )
assigned_clusters <- hm_dat[rownames(hm_dat), "cluster100"]

row_annot <- rowAnnotation(
  Infection = hm_dat$Stimulation,
  # Coculture = hm_dat$Treatment,
  Hour = hm_dat$Timepoint,
  col = list(
    Infection = c("None" = "#44AA99", "HRV" = "#DDCC77"),
    # Coculture = c("Alone" = "#88CCEE", "AEC coculture" = "#AA4499"),
    Hour = c("0" = "grey90", "4" = "grey70", "8" = "grey50",
             "28" = "grey30", "52" = "grey10")
  )
)

hm <-
Heatmap(hm_mat,
        name = "Mean Pearson R",
        # col = colorRamp2(c(0.35,0.45,0.55), 
        #                  c("white", "pink", "red")), #spearman
        col = colorRamp2(c(0.3,0.4,0.5), 
                         c("white", "pink", "red")), #pearson
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        # row_names_side = "left",
        column_names_side = "bottom",
        left_annotation = row_annot, show_row_names = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Get the assigned cluster for this row
          target_cluster <- assigned_clusters[i]
          if (colnames(hm_mat)[j] == target_cluster) {
            grid.text("*", x, y,
                      gp = gpar(col = "black", fontsize = 12, 
                                fontface = "bold"))
          }
        },
        split = hm_dat$Treatment)

png("figures/review_response/Dwyer_pearson_mast_types.png",
    width = 4, height = 6, units = "in", res=300)
draw(hm)
dev.off()

#### dotplot ####
p2 <- dat_markers_mast100 %>% 
  rownames_to_column("gene") %>% 
  left_join(as.data.frame(dat_markers_mc100) %>% 
              rownames_to_column("gene")) %>% 
  pivot_longer(MCTC_1:Transitional, names_to = "cluster") %>% 
  pivot_longer(contains("lib"), names_to = "libID", values_to = "value2") %>% 
  
  ggplot(aes(x=value,y=value2)) +
  geom_point(alpha=0.2) +
  facet_grid(cluster~libID) +
  theme_bw() +
  stat_cor(method="pearson", color="red") +
  stat_cor(method="spearman", color="blue", label.y.npc = "bottom") +
  geom_abline(intercept = 0, slope=1, color="red")
# p2
ggsave("figures/review_response/sample_corr.png", p2,
       width=40, height=10)
