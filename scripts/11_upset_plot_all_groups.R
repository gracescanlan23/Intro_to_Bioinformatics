if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

library(UpSetR)
library(here)

# ----------------------------
# Create output folder
# ----------------------------

output_dir <- here("output", "upset_plot_all_groups")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ----------------------------
# Read DEG result files
# ----------------------------

villous <- read.csv(
  here("output", "villous_vs_normal", "DEG_results_villous_vs_normal.csv"),
  stringsAsFactors = FALSE
)

tubular <- read.csv(
  here("output", "tubular_vs_normal", "DEG_results_tubular_vs_normal.csv"),
  stringsAsFactors = FALSE
)

early_stage <- read.csv(
  here("output", "early_stage_vs_normal", "DEG_results_early_stage_vs_normal.csv"),
  stringsAsFactors = FALSE
)

advanced <- read.csv(
  here("output", "advanced_vs_normal", "DEG_results_advanced_vs_normal.csv"),
  stringsAsFactors = FALSE
)

# ----------------------------
# Filter significant DEGs
# Use the same rule for all groups
# ----------------------------

villous_genes <- unique(
  villous$Gene_ID[
    !is.na(villous$adjusted_p_value) &
    !is.na(villous$log2FC) &
    villous$adjusted_p_value < 0.05 &
    abs(villous$log2FC) > 1
  ]
)

tubular_genes <- unique(
  tubular$Gene_ID[
    !is.na(tubular$adjusted_p_value) &
    !is.na(tubular$log2FC) &
    tubular$adjusted_p_value < 0.05 &
    abs(tubular$log2FC) > 1
  ]
)

early_stage_genes <- unique(
  early_stage$Gene_ID[
    !is.na(early_stage$adjusted_p_value) &
    !is.na(early_stage$log2FC) &
    early_stage$adjusted_p_value < 0.05 &
    abs(early_stage$log2FC) > 1
  ]
)

advanced_genes <- unique(
  advanced$Gene_ID[
    !is.na(advanced$adjusted_p_value) &
    !is.na(advanced$log2FC) &
    advanced$adjusted_p_value < 0.05 &
    abs(advanced$log2FC) > 1
  ]
)

# ----------------------------
# Print DEG counts
# ----------------------------

cat("Villous vs Normal DEGs:", length(villous_genes), "\n")
cat("Tubular vs Normal DEGs:", length(tubular_genes), "\n")
cat("Early Stage vs Normal DEGs:", length(early_stage_genes), "\n")
cat("Advanced vs Normal DEGs:", length(advanced_genes), "\n")

# ----------------------------
# Build list for UpSet plot
# ----------------------------

deg_list <- list(
  Villous = villous_genes,
  Tubular = tubular_genes,
  EarlyStage = early_stage_genes,
  Advanced = advanced_genes
)

# Convert to UpSet input format
upset_data <- fromList(deg_list)

# Save the binary membership table too
write.csv(
  upset_data,
  file.path(output_dir, "upset_membership_table.csv"),
  row.names = TRUE
)

# ----------------------------
# Make and save UpSet plot
# ----------------------------

pdf(file.path(output_dir, "upset_plot_all_groups.pdf"), width = 10, height = 6)

print(
  upset(
    upset_data,
    nsets = 5,
    nintersects = 20,
    order.by = "freq"
  )
)

dev.off()

upset(
  upset_data,
  nsets = 4,
  nintersects = NA,
  order.by = "freq",
  mainbar.y.label = "Shared DEG Intersections",
  sets.x.label = "DEGs per Comparison",
  text.scale = 1.3
)

dev.off()

cat("UpSet plot saved to:", file.path(output_dir, "upset_plot_all_groups.pdf"), "\n")
cat("Membership table saved to:", file.path(output_dir, "upset_membership_table.csv"), "\n")