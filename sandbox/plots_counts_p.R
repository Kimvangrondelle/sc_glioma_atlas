data <- pseudo_diaz1002_with_tau

#sums <- rowSums(diaz_10022_counts)
means <- rowMeans(diaz_10022_counts)

medians <- apply(data, 1, median)
data$max <- apply(data, 1, max)

data$median <- medians
# data$sum <- sums
data$mean <- means


data 
ordered_data <- data %>%
  arrange(rownames(data))
ordered_data$gene <- rownames(ordered_data)


library(ggplot2)
lation <- cor(data$median, data$tau_score_diaz, use = "complete.obs")
print(lation)

ggplot(data, aes(y = median, x = tau_score_diaz)) +
  geom_point(color = "blue", size = 1) +  # Points with color and size
  labs(title = paste("Median counts vs tau score", lation),
       x = "tau score",
       y = "sum of counts") +
  # ylim(2, 1000) +
  theme_minimal()


markers_negbinom_diaz10022 <- FindAllMarkers(pseudo_diaz10022, test.use = "negbinom", min.cells.group = 0, min.cells.feature = 0, 
                                           #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                           return.thresh = 1.1, logfc.threshold = 0, min.pct = 0)
library(dplyr)
markers_negbinom_diaz10022 %>%  dplyr::group_by(cluster) %>% dplyr::select(gene) %>% dplyr::tally()
markers_deseq2_diaz10022 <- FindAllMarkers(pseudo_diaz10022, test.use = "DESeq2", min.cells.group = 0, min.cells.feature = 0, 
                                             #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                             return.thresh = 1.1, logfc.threshold = 0, min.pct = 0)
median_p_val_des <- markers_deseq2_diaz10022 %>%
  group_by(gene) %>%
  summarise(median_p_deseq = median(p_val))
median_p_val_des

min_p_val_des <- markers_deseq2_diaz10022 %>%
  group_by(gene) %>%
  summarise(min_p_deseq = min(p_val))
min_p_val_des

median_p_val_neg <- markers_negbinom_diaz10022 %>%
  group_by(gene) %>%
  summarise(median_p_neg = median(p_val))
median_p_val_neg

min_p_val_neg <- markers_negbinom_diaz10022 %>%
  group_by(gene) %>%
  summarise(min_p_neg = min(p_val))
min_p_val_neg

ggplot(data, aes(y = median, x = tau_score_diaz)) +
  geom_point(color = "blue", size = 1) +  # Points with color and size
  labs(title = paste("Median counts vs tau score", lation),
       x = "tau score",
       y = "sum of counts") +
  # ylim(2, 1000) +
  theme_minimal()

merged_matrix <- merge(median_p_val_neg, ordered_data, by = "gene")
merged_matrix

lation2 <- cor(merged_matrix$median, merged_matrix$median_p_neg, use = "complete.obs")
print(lation2)

ggplot(merged_matrix, aes(y = median, x = median_p_neg)) +
  geom_point(color = "blue", size = 1) +  # Points with color and size
  labs(title = paste("Median counts vs tau score", lation2),
       x = "tau score",
       y = "sum of counts") +
  ylim(0, 10000) +
  theme_minimal()

print(max(merged_matrix))
max_row <- merged_matrix[which.max(merged_matrix$median), ]
max_row

hist(diaz_10022_specificity$tau_score)
indices <- which(diaz_10022_specificity$tau_score < 0.5)
tau_less_half <- diaz_10022_specificity[indices, ]
tau_less_half


min_p_deseq <- pseudo_markers_diaz10022 %>%
  group_by(gene) %>%
  summarise(min_p_deseq = min(-log10(p_val)))
min_p_deseq

diaz_10022_vst$gene <- rownames(diaz_10022_vst)

df2 <- diaz_10022_vst %>% left_join( min_p_deseq, 
                                     by="gene")
df2

ggplot(df2, aes(x=tau_score_vst, y=min_p_deseq)) + 
  geom_point()