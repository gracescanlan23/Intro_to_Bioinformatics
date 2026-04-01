# install if needed
install.packages("patchwork")

library(patchwork)

# Create subfolder for multi-panel figures
if (!dir.exists("output/multi_panel")) {
  dir.create("output/multi_panel")
}
# Volcano Multi-Panel Figure

volcano_plot_advanced <- volcano_plot_advanced + ggtitle("Advanced vs Normal")
volcano_plot_early_stage <- volcano_plot_early_stage + ggtitle("Early-Stage vs Normal")
volcano_plot_tubular <- volcano_plot_tubular + ggtitle("Tubular Adenoma vs Normal")
volcano_plot_villous <- volcano_plot_villous + ggtitle("Villous Adenoma vs Normal")

# Combine volcano plots into 2x2 panel
volcano_panel <- (volcano_plot_advanced | volcano_plot_early_stage) /
                 (volcano_plot_tubular | volcano_plot_villous)

# Add panel labels (A, B, C, D)
volcano_panel <- volcano_panel +
  plot_annotation(tag_levels = "A")

# Display
volcano_panel

# Save image
ggsave("output/multi_panel/volcano_panel.pdf",
       volcano_panel,
       width = 10,
       height = 8,
       dpi = 300)

# MA plot Multi-Panel Figure

ma_plot_advanced <- ma_plot_advanced + ggtitle("Advanced vs Normal")
ma_plot_early_stage <- ma_plot_early_stage + ggtitle("Early-Stage vs Normal")
ma_plot_tubular <- ma_plot_tubular + ggtitle("Tubular Adenoma vs Normal")
ma_plot_villous <- ma_plot_villous + ggtitle("Villous Adenoma vs Normal")

# Combine MA plots into 2x2 panel
ma_panel <- (ma_plot_advanced | ma_plot_early_stage) /
            (ma_plot_tubular | ma_plot_villous)

# Add panel labels (A, B, C, D)
ma_panel <- ma_panel +
  plot_annotation(tag_levels = "A")

# Display
ma_panel

# Save  image
ggsave("output/multi_panel/ma_panel.pdf",
       ma_panel,
       width = 10,
       height = 8,
       dpi = 300)