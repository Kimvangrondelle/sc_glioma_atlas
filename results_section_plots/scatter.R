#file to create scatterplot of bayes score vs zpex score

#load in both data frames
own <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv", row.names = 1)
bayes <- read.csv("output/celltypes/combined_all/spec_bayes.csv", row.names = 1)
bayes <- as.data.frame(-log10(as.matrix(bayes)))

#make both dataframes longer so that each gene appears as times as number of cell types present. -- to make scatterplot and color on cell type
long_own <- own %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = c("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_Astro", "sum_Neuron", "sum_OPC", "sum_Tcell", "sum_ATumor", "sum_ADiv_tumor", "sum_OTumor", "sum_O_notsure"), 
               names_to = "CellType", 
               values_to = "SpecificityScoreOwn")

long_bayes <- bayes %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = c("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_Astro", "sum_Neuron", "sum_OPC", "sum_Tcell", "sum_ATumor", "sum_ADiv_tumor", "sum_OTumor", "sum_O_notsure"), 
               names_to = "CellType", 
               values_to = "SpecificityScoreBayes")

#merge dfs on gene and cell type
merged_df <- inner_join(long_own, long_bayes, by = c("gene", "CellType"))

#plot scores 
ggplot(merged_df, aes(x = SpecificityScoreBayes, y = log(1+SpecificityScoreOwn), color = CellType)) +
  geom_point(pch=19, alpha=0.5) +
  # geom_smooth(method = "loess", se = FALSE) +
  ggpubr::stat_cor(method= "spearman")+ # to get correlation per cell type
  labs(title = "Scatterplot of Own vs. Bayes Specificity Score by Cell Type for all samples combined",
       x = "-log10(Specificity Score Bayes)",
       y = "log(1+Specificity Score Own)") +
  theme_minimal()

# calculate overall spearman correlation between bayes and zpex
own_scores <- log(1+merged_df$SpecificityScoreOwn)
bayes_scores <- merged_df$SpecificityScoreBayes

combined <- as.data.frame(cbind(own_scores, bayes_scores))
combined <- combined %>% replace(is.na(.), 0)

cor(combined$own_scores, combined$bayes_scores, method = "spearman")