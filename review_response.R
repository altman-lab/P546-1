library(tidyverse)
library(limma)
library(ggpubr)
library(patchwork)
library(multcompView)

modsOI <- c(1,2,3,4,9,21)

#### Data ####
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
  filter(Treatment=="Alone") %>% 
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

#### Stats ####
lm_result <- readRDS("dat_clean/WGCNA_output/mc_module_lm_output.rds")
cOI <- data.frame(
  contrast_ref = c("Baseline Alone","Baseline Alone","Baseline Alone","Baseline Alone",
                   "4 hrs Alone","4 hrs Alone","4 hrs Alone",
                   "8 hrs Alone","8 hrs Alone",
                   "28 hrs Alone"),
  contrast_lvl = c("4 hrs Alone","8 hrs Alone","28 hrs Alone","52 hrs Alone",
                   "8 hrs Alone","28 hrs Alone","52 hrs Alone",
                   "28 hrs Alone","52 hrs Alone",
                   "52 hrs Alone")
)

dat_max <- viral_mod_full_df %>% 
  group_by(module) %>% 
  summarise(maxE=max(value)) %>% 
  ungroup()

ct_lvls <- c("Baseline Alone","4 hrs Alone","8 hrs Alone","28 hrs Alone",
             "52 hrs Alone")

lm_result_format <- lm_result$lm.contrast %>% 
  #filter contrasts
  inner_join(cOI) %>% 
  #Fix module names
  mutate(module=gsub("^0","",gene),
         module=as.numeric(module)) %>% 
  select(module, contrast_ref, contrast_lvl, FDR) %>%
  #order contrasts
  # mutate(contrast_ref = factor(contrast_ref, levels=ct_lvls),
  #        contrast_lvl = factor(contrast_lvl, levels=ct_lvls)) %>% 
  mutate(contrast_ref = gsub(" Alone","",contrast_ref),
         contrast_lvl = gsub(" Alone","",contrast_lvl)) %>% 
  
  rename(group1=contrast_ref,group2=contrast_lvl) %>% 
  #format FDR labels
  mutate(FDR = signif(FDR, digits=2),
         `.y.`= "value",
         signif = case_when(FDR<0.01~"***",
                            FDR<0.05~"**",
                            FDR<0.2~"*")) %>% 
  #y position
  left_join(dat_max) %>% 
  arrange(module, group1, group2) %>% 
  group_by(module) %>% 
  mutate(rank=row_number(),
         y.position=maxE+rank*0.3) %>% 
  select(`.y.`, group1, group2, FDR, signif, y.position, module)


#### plot ####
plot.ls <- list() #FDR bars
plot.ls2 <- list() #signif letters
for(m in modsOI){
  print(m)
  dat_lines_temp <- dat_lines %>% 
    filter(module==m)
  df_temp <- viral_mod_full_df %>% 
    filter(module == m)
  stat_temp <- lm_result_format %>% 
    #Signif cutoff
    filter(FDR<0.2) %>% 
    filter(module == m)
  
  base_plot <- ggplot(df_temp, aes(x = Timepoint, y = value)) +
    geom_boxplot(aes(color = Treatment)) +
    geom_point(position = position_dodge(width = 0.75), 
               aes(group = Treatment), size = 1) +
    geom_line(data = dat_lines_temp,
              mapping = aes(x = Timepoint, y = avg, color=Treatment, group=1),
              linetype = "solid", size = 1) + 
    stat_summary(fun = mean, geom = "point", 
                 aes(shape = Treatment, fill = Treatment), 
                 size = 3.00, position = "identity") +
    scale_shape_manual(values = c(21,24)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black")) +
    labs(y=paste("Mean Module",m,"Expression (log2 CPM)"),
         subtitle = paste("Module",m)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = "none")
  
  plot.ls[[as.character(m)]] <- base_plot +
    stat_pvalue_manual(data = stat_temp, label="FDR")
  
  #Get significance letters
  stat_temp2 <- lm_result_format %>% 
    filter(module == m)
  signif_vec <- stat_temp2$FDR
  names(signif_vec) <- paste(stat_temp2$group1,
                             stat_temp2$group2,sep="-")
  signif_letters <- multcompLetters(
    x = signif_vec,
    compare = "<",
    threshold = 0.2,
    Letters = c(letters, LETTERS, "."),
    reversed = FALSE
  )$Letters %>% 
    as.data.frame() %>% 
    rownames_to_column("Timepoint") %>% 
    mutate(y=max(df_temp$value+0.2))
  
  plot.ls2[[as.character(m)]] <- base_plot +
    geom_text(data = signif_letters,
              aes(x = Timepoint, y = y, label = `.`),
              inherit.aes = FALSE,
              vjust = 0) 
}

p <- wrap_plots(plot.ls)
p2 <- wrap_plots(plot.ls2)
# p
ggsave("figures/review_response/modules_mc_alone_time_v2.png", p,
       width=12, height=8)
ggsave("figures/review_response/modules_mc_alone_time_alt.png", p2,
       width=12, height=8)
