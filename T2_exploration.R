library(tidyverse)
library(limma)
library(RNAetc)
# T2 gene set (average of IL4, IL5, and IL13) 
# see if this correlates with Module expression 
#pull the gene level data from each of the genes (IL4, IL5, IL13, IL1B, IL33, IFNB1)

#### Define genes of interest ####

geneOI <- c("IL4", "IL5", "IL13", "IL1B", "IL33", "IFNB1")
load("dat_clean/ensembl_gene_dictionary.RData")
ens_gene_dict <- ens_gene_dict %>% 
  filter(geneName %in% geneOI)

#### AEC data ####
load("dat_clean/voom_objects/aec_voom.RData")
aec_voom$genes <- NULL

dat_aec <- collapse_voom(aec_voom, libraryID = "libID",
                         geneID = "ensembl_gene_id")
dat_aec <- dat_aec %>% 
  inner_join(ens_gene_dict, 
             by = join_by(ensembl_gene_id)) %>% 
  select(ensembl_gene_id, geneName, libID, value, 
         Study.Group, Treatment, Stimulation, Timepoint)

unique(dat_aec$geneName)
# "IL1B" "IL33"

#### MC data ####
load("dat_clean/voom_objects/mc_voom_viral.RData")
dat_mc <- collapse_voom(mc_voom_viral, libraryID = "libID",
                         geneID = "ensembl_gene_id")
dat_mc <- dat_mc %>% 
  inner_join(ens_gene_dict,
             by = join_by(ensembl_gene_id)) %>% 
  select(ensembl_gene_id, geneName, libID, value, 
         Study.Group, Treatment, Stimulation, Timepoint, ptID)

unique(dat_mc$geneName)
# "IL5"   "IL1B"  "IL33"  "IL13"  "IFNB1"

#### Rare genes ####
aec_rare <- read_csv("dat_clean/rare_aec_genes.csv") %>% 
  filter(ensembl_gene_id %in% ens_gene_dict$ensembl_gene_id)
mc_rare <- read_csv("dat_clean/rare_mc_genes.csv") %>% 
  filter(ensembl_gene_id %in% ens_gene_dict$ensembl_gene_id)

rare_all <- aec_rare %>% 
  mutate(cell="AEC",.before = 0) %>% 
  bind_rows(mc_rare %>% 
              mutate(cell="MC",.before = 0)) 
rare_all %>% 
  write_csv("~/Desktop/P546-1_rare_genes.csv")

#### Raw counts ####
raw_aec <- read_csv("dat_raw/P546-1_AEC_raw_counts.csv") %>% 
  column_to_rownames("ensembl_gene_id")
cpm_aec <- edgeR::cpm(raw_aec, log=FALSE) %>% 
  as.data.frame() %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  filter(ensembl_gene_id %in% ens_gene_dict$ensembl_gene_id)

raw_mc <- read_csv("dat_raw/P546-1_MC_raw_counts.csv") %>% 
  column_to_rownames("ensembl_gene_id")
cpm_mc <- edgeR::cpm(raw_mc, log=FALSE) %>% 
  as.data.frame() %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  filter(ensembl_gene_id %in% ens_gene_dict$ensembl_gene_id)

cpm_all <- cpm_aec %>% 
  pivot_longer(-ensembl_gene_id, names_to = "libID", values_to = "cpm") %>% 
  left_join(aec_voom$targets %>%
              select(libID, Study.Group, Treatment, Stimulation, Timepoint),
            by = join_by(libID)) %>% 
  bind_rows(
    cpm_mc %>% 
      pivot_longer(-ensembl_gene_id, names_to = "libID", values_to = "cpm") %>% 
      left_join(mc_voom_viral$targets %>%
                  select(libID, Study.Group, Treatment, Stimulation, Timepoint),
                by = join_by(libID))
  ) %>% 
  left_join(ens_gene_dict, by = join_by(ensembl_gene_id)) 

p1 <- cpm_all %>% 
  mutate(Stimulation = fct_recode(Stimulation, "HRV"="HRV Infection",
                                  "HRV MCs"="HRV-Infected MCs")) %>% 
  mutate(x=paste(Stimulation, Timepoint, sep="\n"),
         Treatment = case_when(Study.Group=="AEC" & Treatment=="Alone"~"AEC alone",
                               Study.Group=="AEC" & Treatment=="MCs"~"AEC with MC",
                               Study.Group=="MC" & Treatment=="Alone"~"MC with AEC",
                               Study.Group=="MC" & Treatment=="AECs"~"MC alone")) %>% 
  left_join(rare_all %>% rename(Study.Group=cell) %>% 
              select(Study.Group, ensembl_gene_id) %>% 
              mutate(col.group="FAIL")) %>% 
  mutate(col.group=ifelse(is.na(col.group), "PASS",col.group)) %>% 
  ggplot(aes(x=x, y=cpm, color=col.group)) +
  geom_jitter(height=0,width=0.2) +
  theme_bw(base_size = 9) +
  facet_wrap(Treatment~geneName, scales="free", ncol=6) +
  labs(y="CPM",x="", color="Gene filter") +
  scale_color_manual(values=c("FAIL"="red","PASS"="black"))

# p1
ggsave("~/Desktop/select_gene_CPM.png", p1, width=12, height=8)
