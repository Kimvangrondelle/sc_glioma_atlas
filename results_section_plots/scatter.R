own <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv", row.names = 1)


bayes <- read.csv("output/celltypes/combined_all/spec_bayes.csv", row.names = 1)
bayes <- as.data.frame(-log10(as.matrix(bayes)))

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

merged_df <- inner_join(long_own, long_bayes, by = c("gene", "CellType"))

ggplot(merged_df, aes(x = SpecificityScoreBayes, y = log(1+SpecificityScoreOwn), color = CellType)) +
  geom_point(pch=19, alpha=0.5) +
  # geom_smooth(method = "loess", se = FALSE) +
  ggpubr::stat_cor(method= "spearman")+
  labs(title = "Scatterplot of Own vs. Bayes Specificity Score by Cell Type for all samples combined",
       x = "-log10(Specificity Score Bayes)",
       y = "log(1+Specificity Score Own)") +
  theme_minimal()

own_scores <- log(1+merged_df$SpecificityScoreOwn)
bayes_scores <- merged_df$SpecificityScoreBayes

combined <- as.data.frame(cbind(own_scores, bayes_scores))
combined <- combined %>% replace(is.na(.), 0)

cor(combined$own_scores, combined$bayes_scores, method = "spearman")