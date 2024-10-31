nebula_df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_nebula_gsea_results.csv")
deseq2_df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_deseq2_gsea_results.csv")
esvd_df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_esvd_gsea_results.csv")

rownames(nebula_df) <- nebula_df$ID
nebula_df <- nebula_df[esvd_df$ID,] 

rownames(deseq2_df) <- deseq2_df$ID
deseq2_df <- deseq2_df[esvd_df$ID,] 

logpvalue_esvd <- -log10(esvd_df$pvalue)
names(logpvalue_esvd) <- esvd_df$ID
logpvalue_nebula <- -log10(nebula_df$pvalue)
names(logpvalue_nebula) <- nebula_df$ID
logpvalue_deseq2 <- -log10(deseq2_df$pvalue)
names(logpvalue_deseq2) <- deseq2_df$ID

df <- data.frame(id = esvd_df$ID,
                 Description = esvd_df$Description,
                 esvd = logpvalue_esvd,
                 nebula = logpvalue_nebula,
                 deseq2 = logpvalue_deseq2)
combined_df <- df

library(tidyr)

df_long <- df %>%
  pivot_longer(
    cols = c(esvd, nebula, deseq2),
    names_to = "method",
    values_to = "logpvalue"
  )

library(tidyverse)
write_csv(combined_df, 
          "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_combined_df_gsea.csv")

write_csv(df_long, 
          "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_df_long_gsea.csv")

########################################################################
######################## plot horizontal barplot ##########################################
########################################################################

df_long$method <- factor(df_long$method, levels = c("deseq2", "nebula", "esvd"))
plot1 <- ggplot(df_long, aes(x = Description, y = logpvalue, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" to compare side-by-side bars
  labs(x = "Pathways", y = "-log10(p-value)", 
       title = "-log10(p-value) of Gene Pathways") + scale_fill_manual(values = c(
         "deseq2" = "#A0D3FF",  # Light blue
         "nebula" = "#FF8AD0",  # Pink
         "esvd" = "#A886F6"     # Purple
       ))+
  theme_minimal() 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup8_gsea_horizontal_barplot.png"),
                height = 900, width = 3000, units = "px")
#######
plot2 <- ggplot(df_long, aes(x = Description, y = logpvalue, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" to compare side-by-side bars
  labs(x = "Pathways", y = "-log10(p-value)", 
       title = "-log10(p-value) of Gene Pathways") + scale_fill_manual(values = c(
         "deseq2" = "#A0D3FF",  # Light blue
         "nebula" = "#FF8AD0",  # Pink
         "esvd" = "#A886F6"     # Purple
       ))+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(stringr)
plot2 <- plot2 + scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 5
    ),
    plot.margin = margin(b = 40, l = 10, r = 10, t = 10)
  )

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/"
ggplot2::ggsave(plot2, 
                filename = paste0(plot_folder, "Writeup8_gsea_horizontal_barplot_test.png"),
                height = 900, width = 3000, units = "px")

