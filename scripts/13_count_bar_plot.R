library(ggplot2)
library(tidyr)
library(dplyr)

# Create DEG counts dataframe
deg_counts <- data.frame(
  Tumor = c("Advanced Tumor", "Early-Stage Tumor", "Tubular Adenoma", "Villous Adenoma"),
  Up = c(advanced_up_count, early_up_count, tubular_up_count, villous_up_count),
  Down = c(advanced_down_count, early_down_count, tubular_down_count, villous_down_count)
)

# Convert to long format
deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c("Up", "Down"),
               names_to = "Regulation",
               values_to = "Count")

# Total counts for labels
deg_totals <- deg_counts_long %>%
  group_by(Tumor) %>%
  summarise(Total = sum(Count), .groups = "drop")

# Plot
deg_count_barplot <- ggplot(deg_counts_long, aes(x = Tumor, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  
  # Total labels on top
  geom_text(data = deg_totals,
            aes(x = Tumor, y = Total, label = Total),
            inherit.aes = FALSE,
            vjust = -0.4,
            size = 4) +
  
  theme_minimal() +
  labs(
    title = "Differentially Expressed Genes Across Tumor Types",
    x = "Tumor Type",
    y = "Number of DEGs"
  ) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display
deg_count_barplot

ggsave(
  filename = "output/count_bar_plot/deg_count_barplot.png",
  plot = deg_count_barplot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "output/count_bar_plot/deg_count_barplot.pdf",
  plot = deg_count_barplot,
  width = 8,
  height = 6
)