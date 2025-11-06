library(tidyverse)
library(limma)
library(patchwork)
modsOI <- c(1,2,3,4,9,21)

#### Expression data ####
#Mean module expression
mean_mod_E <- read.csv(file = "dat_clean/WGCNA_output/mc_mean_module_exprs.csv") %>% 
  dplyr::select(-X)
#voom object
load("dat_clean/voom_objects/mc_voom.RData")

#combine
viral_mod_full_df <- mean_mod_E  %>% 
  #modules of interest
  filter(module %in% modsOI) %>% 
  #transpose
  pivot_longer(-module, names_to = "libID") %>% 
  #add metadata
  left_join(mc_voom$targets, by = "libID") %>% 
  #just mast cells alone
  # filter(Treatment=="Alone") %>% 
  #Create concat variable
  mutate(TimeStimTreat = paste0(Timepoint, "_", Stimulation, "_", Treatment), 
         .after=Timepoint) %>%
  mutate(TimeStimTreat = factor(TimeStimTreat, 
                                levels = c("Baseline_None_Alone",
                                           "4 hrs_HRV Infection_Alone",
                                           "8 hrs_HRV Infection_Alone",
                                           "8 hrs_HRV Infection_AECs",
                                           "28 hrs_HRV Infection_Alone",
                                           "28 hrs_HRV Infection_AECs",
                                           "52 hrs_HRV Infection_Alone",
                                           "52 hrs_HRV Infection_AECs")))

#Mean trend lines
dat_lines <- viral_mod_full_df %>% 
  group_by(module, Timepoint, Treatment, TimeStimTreat) %>% 
  summarise(avg=mean(value)) %>% 
  ungroup()

#### Boxplots ####
plot.ls <- list() 
for(m in modsOI){
  print(m)
  dat_lines_temp <- dat_lines %>% 
    filter(module==m)
  df_temp <- viral_mod_full_df %>% 
    filter(module == m)
  ylims <- c(min(df_temp$value), max(df_temp$value)+0.2)
  
  base_plot <- ggplot(df_temp, aes(x = Timepoint, y = value)) +
    geom_boxplot(aes(color = Treatment), show.legend = FALSE, 
                 outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.75), 
               aes(group = Treatment), size = 1) +
    geom_line(data = dat_lines_temp %>% 
                filter(TimeStimTreat %in% c("Baseline_None_Alone",
                                            "4 hrs_HRV Infection_Alone")),
              mapping = aes(x = Timepoint, y = avg, color=Treatment, group=1),
              color = "blue", linetype = "solid", size = 1) +
    geom_line(data = dat_lines_temp %>% 
                filter(TimeStimTreat %in% c("4 hrs_HRV Infection_Alone",
                                            "8 hrs_HRV Infection_Alone",
                                            "28 hrs_HRV Infection_Alone",
                                            "52 hrs_HRV Infection_Alone")),
              mapping = aes(x = Timepoint, y = avg, color=Treatment, group=1),
              color = "#F8766D", linetype = "solid", size = 1) +
    geom_line(data = dat_lines_temp %>% 
                filter(TimeStimTreat %in% c("4 hrs_HRV Infection_Alone",
                                            "8 hrs_HRV Infection_AECs",
                                            "28 hrs_HRV Infection_AECs",
                                            "52 hrs_HRV Infection_AECs")),
              mapping = aes(x = Timepoint, y = avg, color=Treatment, group=1),
              color = "#00BFC4", linetype = "solid", size = 1) +
    stat_summary(fun = mean, geom = "point",
                 aes(shape = Treatment, fill = Treatment),
                 size = 2, position = "identity") +
    scale_shape_manual(values = c(21,24)) +
    labs(y=paste("Mean Expression (log2 CPM)"),
         subtitle = paste("Module",m),
         fill="Co-culture",shape="Co-culture") +
    theme_bw(base_size = 9) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black"),
          plot.subtitle = element_text(hjust = 0.5))  +
    guides(
      fill = guide_legend(nrow=2),
      shape = guide_legend(nrow=2)
    ) +
    lims(y=ylims)
  
  if(m=="1"){
    plot.ls[[as.character(m)]] <- base_plot+ theme(legend.position = "bottom")
  }else{
    plot.ls[[as.character(m)]] <- base_plot + theme(legend.position = "none")
  }
}
#
#### Enrich data ####
enr_mod_GO <- readxl::read_excel(path = "dat_clean/WGCNA_output/module_enrichment_GOBP.xlsx") %>% 
  filter(group %in% modsOI) %>% 
  mutate(pathway = gsub("GOBP_", "", pathway)) %>% 
  mutate(pathway = tolower(gsub("_", " ", pathway))) %>% 
  filter(FDR<0.1)

#Select pathways
enr_key <- data.frame(
  group = "1",
  pathway = c("response to virus", 
              "regulation of innate immune response",
              "viral genome replication", 
              "response to interferon gamma",
              "response to type i interferon",
              "type i interferon production"),
  pw = c("response to virus",
         "regulation of innate immune response",
         "viral genome replication", 
         "response to IFNG",
         "response to Type I IFN",
         "Type I IFN production")) %>% 
  bind_rows(data.frame(
    group = "2",
    pathway = c("dna replication", 
                "chromosome segregation", 
                "mitotic nuclear division",
                "mitotic sister chromatid segregation"),
    pw = c("DNA replication", 
           "chromosome segregation", 
           "mitotic nuclear division",
           "mitotic sister chromatid segregation"))) %>% 
  bind_rows(data.frame(
    group = "3",
    pathway = c("mast cell activation",
                "regulation of mast cell activation",
                "leukocyte degranulation",
                "regulation of leukocyte degranulation",
                "myeloid leukocyte mediated immunity",
                "myeloid leukocyte activation"),
    pw = c("mast cell activation",
           "regulation of mast cell activation",
           "leukocyte degranulation",
           "regulation of leukocyte degranulation",
           "myeloid leukocyte mediated immunity",
           "myeloid leukocyte activation"))) %>% 
  bind_rows(data.frame(
    group = "4",
    pathway = c("regulation of apoptotic signaling pathway",
                #"extrinsic apoptotic signaling pathway",
                #"intrinsic apoptotic signaling pathway",
                "apoptotic mitochondrial changes",
                "transforming growth factor beta receptor signaling pathway",
                "response to transforming growth factor beta"),
    pw = c("regulation of apoptotic signaling pathway",
           #"extrinsic apoptotic signaling pathway",
           #"intrinsic apoptotic signaling pathway",
           "apoptotic mitochondrial changes",
           "TGF-beta receptor signaling pathway",
           "response to TGF-beta"))) %>% 
  bind_rows(data.frame(
    group = "9",
    pathway = c("protein folding",
                "atp biosynthetic process", 
                "ribosome assembly", 
                "regulation of cellular protein catabolic process"),
    pw = c("protein folding",
           "ATP biosynthetic process", 
           "ribosome assembly", 
           "regulation of cellular protein catabolic process"))) %>% 
  bind_rows(data.frame(
    group = "21",
    pathway = c("atp metabolic process", 
                "atp synthesis coupled electron transport", 
                "aerobic respiration",
                "cellular respiration", 
                "nadh dehydrogenase complex assembly"),
    pw = c("ATP metabolic process", 
           "ATP synthesis coupled electron transport", 
           "aerobic respiration",
           "cellular respiration", 
           "NADH dehydrogenase complex assembly")))

#clean names fxn
add_hard_returns <- function(text, width = 80) {
  wrapped <- stringi::stri_wrap(text, width = width, simplify = FALSE, use_length = TRUE)
  vapply(wrapped, function(x) paste(x, collapse = "\n"), character(1))
}

#Filter data and clean
enr_selected <- enr_mod_GO %>% 
  inner_join(enr_key) %>%
  mutate(`FDR Bracket` = case_when(FDR < 0.001 ~ "FDR < 0.001",
                                   FDR < 0.01 ~ "FDR < 0.01",
                                   FDR < 0.1 ~ "FDR < 0.1")) %>%
  mutate(`FDR Bracket` = factor(`FDR Bracket`,
                                levels = c("FDR < 0.001", 
                                           "FDR < 0.01",
                                           "FDR < 0.1"))) %>% 
  mutate(pw_clean =add_hard_returns(pw, width = 23))

#### Enrich plots ####
plot.ls2 <- list() 

for(m in modsOI){
  print(m)
  enr_selected_temp <- enr_selected %>% 
    filter(group == m)
  
  plot2 <-
    ggplot(enr_selected_temp,
           aes(y = fct_reorder2(pw_clean, `FDR Bracket`,-group_in_pathway),
               x = 1,
               size = group_in_pathway,
               color = `FDR Bracket`)) +
    geom_point() +
    # dummy layer: only color mapped
    # geom_point(
    #   data = tibble(`FDR Bracket` = c("FDR < 0.001","FDR < 0.01",
    #                                   "FDR < 0.1")),
    #   aes(color = `FDR Bracket`, x=1,y=1),
    #   inherit.aes = FALSE, size = 0
    # ) +
    scale_color_manual(values = c("FDR < 0.001" = "#4DAF4A",
                                  "FDR < 0.01"  = "#984EA3",
                                  "FDR < 0.1"   = "#F781BF"), 
                       drop = FALSE,
                       breaks=c("FDR < 0.001","FDR < 0.01","FDR < 0.1")) +
    scale_size(range = c(0, 5), breaks = c(5,25,50,100), 
               limits = c(0,100)) +
    theme_bw(base_size = 9) +
    labs(size='Number of Genes\nin Pathway') +
    theme(axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    guides(
      color = guide_legend(override.aes = list(size = 3), nrow=3,
                           keyheight = unit(0.1, "cm")),
      size = guide_legend(nrow=2))
  
  if(m=="3"){
    plot.ls2[[as.character(m)]] <- plot2 + theme(legend.position = "bottom")
  }else{
    plot.ls2[[as.character(m)]] <- plot2 + theme(legend.position = "none")
  }
}


#### Combine ####
lo <- "
AAAAAABCCCCCCD
EEEEEEFGGGGGGH
IIIIIIJKKKKKKL
"
p <- plot.ls[[1]] + plot.ls2[[1]]+
  plot.ls[[2]] + plot.ls2[[2]]+
  plot.ls[[3]] + plot.ls2[[3]]+
  plot.ls[[4]] + plot.ls2[[4]]+
  plot.ls[[5]] + plot.ls2[[5]]+
  plot.ls[[6]] + plot.ls2[[6]]+
  
  plot_annotation(tag_levels = list(c("A","","B","","C","",
                                      "D","","E","","F","")),
                  theme = theme(legend.position = "bottom")) +
  plot_layout(guides = "collect", design = lo)

# p
ggsave("figures/Figure4_all_mods.png", p,
       width=8, height=7)
ggsave("figures/Figure4_all_mods.pdf", p,
       width=8, height=7)

