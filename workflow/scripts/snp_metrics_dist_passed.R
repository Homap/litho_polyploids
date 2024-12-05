#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the distributions of the INFO fields
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(gridExtra)
library(ggplot2)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load info field
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
info_field <- fread("results/filtered_vcf/info_field_passed.txt", header = FALSE, sep = "\t",
                          col.names = c("CHROM","POS","AF","BaseQRankSum","DP","ExcessHet","FS","InbreedingCoeff","MQ","MQRankSum","QD","ReadPosRankSum","SOR"))
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# To get the distributions, I will use a random sample representing 1% of the data
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set seed for reproducibility
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sample 1% of the rows
n_rows <- nrow(info_field)
n_sample <- ceiling(n_rows * 0.05)
# Sample and save to a new data.table
sampled_indices <- sample(1:n_rows, n_sample)   # This gives the row indices of the sampled rows
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get 1% of the info field
sampled_info_field <- info_field[sampled_indices]       # Get the sampled rows
# Convert columns to numeric where necessary
sampled_info_field[] <- lapply(sampled_info_field, function(x) as.numeric(as.character(x)))
# Define a function to create density plots for each field
plot_density <- function(data, column) {
  ggplot(data, aes_string(x = column)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = paste("Distribution of", column), x = column, y = "Density") +
    theme_minimal()
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# List of columns to plot
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
columns_to_plot <- c("AF","BaseQRankSum","DP","ExcessHet","FS","InbreedingCoeff","MQ","MQRankSum","QD","ReadPosRankSum","SOR")
# Set up PDF output
pdf("results/filtered_vcf/info_plots_grid_passed.pdf", width = 11, height = 8.5) # Adjust width and height as needed
# Generate each column's density plot
plot_list <- lapply(columns_to_plot, function(col) plot_density(sampled_info_field, col))
# Arrange plots in a 3-column grid and output to PDF
do.call(grid.arrange, c(plot_list, ncol = 3))
# Close the PDF device
dev.off()