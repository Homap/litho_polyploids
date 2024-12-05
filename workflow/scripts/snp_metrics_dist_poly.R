#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the distributions of the INFO fields
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(gridExtra)
library(ggplot2)
library(reshape2)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sample names
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
diploid_sample_names <- c("P18758_101_S1_L001", "P18758_102_S2_L001", "P18758_103_S3_L001", "P18758_104_S4_L001", 
"P18758_131_S31_L001", "P18758_132_S32_L001", "P18758_133_S33_L002", "P18758_134_S34_L002", 
"P18758_161_S61_L002", "P18758_162_S62_L002", "P18758_163_S63_L002", "P18758_164_S64_L002")
tetraploid_sample_names <- c("P18758_169_S79_L004", "P18758_170_S80_L004", "P18758_171_S81_L004", "P18758_172_S82_L004", "P18758_173_S83_L004", "P18758_174_S84_L004")
hexaploid_sample_names <- c("P18758_175_S85_L004", "P18758_176_S86_L004", "P18758_177_S87_L004", "P18758_178_S88_L004")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read the file with total mapped reads for diploids, tetraploids and hexaploids
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
total_reads_sample_diploid <- read.table("config/total_mapped_reads.txt", header = T)[c(1:4, 31:34, 61:64),]
total_reads_sample_tetraploid <- read.table("config/total_mapped_reads.txt", header = T)[69:74,]
total_reads_sample_hexaploid <- read.table("config/total_mapped_reads.txt", header = T)[75:78,]
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load info field
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
info_field <- fread("results/filtered_vcf_poly/info_field.txt", header = FALSE, sep = "\t",
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
n_sample <- ceiling(n_rows * 0.01)
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
pdf("results/filtered_vcf_poly/info_plots_grid.pdf", width = 11, height = 8.5) # Adjust width and height as needed
# Generate each column's density plot
plot_list <- lapply(columns_to_plot, function(col) plot_density(sampled_info_field, col))
# Arrange plots in a 3-column grid and output to PDF
do.call(grid.arrange, c(plot_list, ncol = 3))
# Close the PDF device
dev.off()
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the distributions for DP and normalised DP for diploids, tetraploids and hexaploids
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load the data
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
diploid <- fread("results/filtered_vcf_poly/depth_diploid.txt", sep = "\t", header = T)
tetraploid <- fread("results/filtered_vcf_poly/depth_tetraploid.txt", sep = "\t", header = T)
hexaploid <- fread("results/filtered_vcf_poly/depth_hexaploid.txt", sep = "\t", header = T)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sample and save to a new data.table
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
sampled_diploid <- diploid[sampled_indices]       # Get the sampled rows
sampled_tetraploid <- tetraploid[sampled_indices]       # Get the sampled rows
sampled_hexploid <- hexaploid[sampled_indices]       # Get the sampled rows

colnames(sampled_diploid) <- c("contig", "position", "total_depth", diploid_sample_names)
colnames(sampled_tetraploid) <- c("contig", "position", "total_depth", tetraploid_sample_names)
colnames(sampled_hexploid) <- c("contig", "position", "total_depth", hexaploid_sample_names)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Write the sampled data
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
fwrite(sampled_diploid, "results/filtered_vcf_poly/sampled_diploid.csv")
fwrite(sampled_tetraploid, "results/filtered_vcf_poly/sampled_tetraploid.csv")
fwrite(sampled_hexploid, "results/filtered_vcf_poly/sampled_hexaploid.csv")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the sampled data
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
sampled_diploid_data <- fread("results/filtered_vcf_poly/sampled_diploid.csv", header=T)
sampled_tetraploid_data <- fread("results/filtered_vcf_poly/sampled_tetraploid.csv", header=T)
sampled_hexaploid_data <- fread("results/filtered_vcf_poly/sampled_hexaploid.csv", header=T)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get a copy of sampled data since I re-write the original later with normalised values
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
sampled_diploid_data_copy <- sampled_diploid_data
sampled_tetraploid_data_copy <- sampled_tetraploid_data
sampled_hexaploid_data_copy <- sampled_hexaploid_data
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Normalize each sample column containing read depth per site by the total reads mapping for that sample
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
for (sample in total_reads_sample_diploid$sample) {
  # Find the total depth for the current sample
  depth_diploid <- total_reads_sample_diploid$totalreads[total_reads_sample_diploid$sample == sample]
  # Divide the corresponding column in df2 by the depth
  sampled_diploid_data[[sample]] <- sampled_diploid_data[[sample]] / depth_diploid
}
for (sample in total_reads_sample_tetraploid$sample) {
  # Find the total depth for the current sample
  depth_tetraploid <- total_reads_sample_tetraploid$totalreads[total_reads_sample_tetraploid$sample == sample]
  # Divide the corresponding column in df2 by the depth
  sampled_tetraploid_data[[sample]] <- sampled_tetraploid_data[[sample]] / depth_tetraploid
}
for (sample in total_reads_sample_hexaploid$sample) {
  # Find the total depth for the current sample
  depth_hexaploid <- total_reads_sample_hexaploid$totalreads[total_reads_sample_hexaploid$sample == sample]
  # Divide the corresponding column in df2 by the depth
  sampled_hexaploid_data[[sample]] <- sampled_hexaploid_data[[sample]] / depth_hexaploid
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the normalized mean of each row, ignoring the first three columns (contig, position, totaldepth)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
diploid_row_means <- rowMeans(sampled_diploid_data[, 4:ncol(sampled_diploid_data)], na.rm = TRUE)
tetraploid_row_means <- rowMeans(sampled_tetraploid_data[, 4:ncol(sampled_tetraploid_data)], na.rm = TRUE)
hexaploid_row_means <- rowMeans(sampled_hexaploid_data[, 4:ncol(sampled_hexaploid_data)], na.rm = TRUE)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# For each position, get a ratio of diploid mean over tetraploid mean
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
diploid_tetraploid_ratio <- diploid_row_means/tetraploid_row_means
diploid_tetraploid_ratio_cleaned <- diploid_tetraploid_ratio[is.finite(diploid_tetraploid_ratio)]
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# For each position, get a ratio of diploid mean over hexaploid mean
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
diploid_hexaploid_ratio <- diploid_row_means/hexaploid_row_means
diploid_hexaploid_ratio_cleaned <- diploid_hexaploid_ratio[is.finite(diploid_hexaploid_ratio)]
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# raw ratio
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
dip_tetra_raw_ratio <- rowMeans(sampled_diploid_data_copy[, 4:ncol(sampled_diploid_data_copy)], na.rm = TRUE)/rowMeans(sampled_tetraploid_data_copy[, 4:ncol(sampled_tetraploid_data_copy)], na.rm = TRUE)
dip_tetra_raw_ratio_cleaned <- dip_tetra_raw_ratio[is.finite(dip_tetra_raw_ratio)]
dip_hexa_raw_ratio <- rowMeans(sampled_diploid_data_copy[, 4:ncol(sampled_diploid_data_copy)], na.rm = TRUE)/rowMeans(sampled_hexaploid_data_copy[, 4:ncol(sampled_hexaploid_data_copy)], na.rm = TRUE)
dip_hexa_raw_ratio_cleaned <- dip_hexa_raw_ratio[is.finite(dip_hexa_raw_ratio)]
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Write a summary table
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate statistics
diploid_stats <- c(mean(rowMeans(sampled_diploid_data_copy[, 4:ncol(sampled_diploid_data_copy)], na.rm = TRUE)), 
                median(rowMeans(sampled_diploid_data_copy[, 4:ncol(sampled_diploid_data_copy)], na.rm = TRUE)), 
                sd(rowMeans(sampled_diploid_data_copy[, 4:ncol(sampled_diploid_data_copy)], na.rm = TRUE)))
tetraploid_stats <- c(mean(rowMeans(sampled_tetraploid_data_copy[, 4:ncol(sampled_tetraploid_data_copy)], na.rm = TRUE)), 
                median(rowMeans(sampled_tetraploid_data_copy[, 4:ncol(sampled_tetraploid_data_copy)], na.rm = TRUE)), 
                sd(rowMeans(sampled_tetraploid_data_copy[, 4:ncol(sampled_tetraploid_data_copy)], na.rm = TRUE)))
hexaploid_stats <- c(mean(rowMeans(sampled_hexaploid_data_copy[, 4:ncol(sampled_hexaploid_data_copy)], na.rm = TRUE)), 
                median(rowMeans(sampled_hexaploid_data_copy[, 4:ncol(sampled_hexaploid_data_copy)], na.rm = TRUE)), 
                sd(rowMeans(sampled_hexaploid_data_copy[, 4:ncol(sampled_hexaploid_data_copy)], na.rm = TRUE)))
# Create a data frame to store the results
stats_table <- data.frame(
  Ploidy = c("Diploid", "Tetraploid", "Hexaploid"),
  Mean = c(diploid_stats[1], tetraploid_stats[1], hexaploid_stats[1]),
  Median = c(diploid_stats[2], tetraploid_stats[2], hexaploid_stats[2]),
  SD = c(diploid_stats[3], tetraploid_stats[3], hexaploid_stats[3])
)
# Write the table to a CSV file
write.csv(stats_table, file = "results/filtered_vcf_poly/summary_statistics_with_ploidy.csv", row.names = FALSE)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set up the PDF output
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf("results/filtered_vcf_poly/diploid_tetra_hexa_ratio.pdf", width = 8.5, height = 11)
# Create density plots for each dataset
p1 <- ggplot(data.frame(value = rowMeans(sampled_diploid_data_copy[, 4:ncol(sampled_diploid_data_copy)], na.rm = TRUE)), aes(x = value)) + 
  geom_density() + 
  ggtitle("Diploid Depth") + 
  theme_minimal()
p2 <- ggplot(data.frame(value = rowMeans(sampled_tetraploid_data_copy[, 4:ncol(sampled_tetraploid_data_copy)], na.rm = TRUE)), aes(x = value)) + 
  geom_density() + 
  ggtitle("Tetraploid Depth") + 
  theme_minimal()
p3 <- ggplot(data.frame(value = rowMeans(sampled_hexaploid_data_copy[, 4:ncol(sampled_hexaploid_data_copy)], na.rm = TRUE)), aes(x = value)) + 
  geom_density() + 
  ggtitle("Hexaploid Depth") + 
  theme_minimal()
p4 <- ggplot(data.frame(value = dip_tetra_raw_ratio_cleaned), aes(x = value)) + 
  geom_density() + 
  ggtitle("Diploid - Tetraploid Ratio") + 
  xlim(0, 2) + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  theme_minimal()
p5 <- ggplot(data.frame(value = dip_hexa_raw_ratio_cleaned), aes(x = value)) + 
  geom_density() + 
  ggtitle("Diploid - Hexaploid Ratio") + 
  xlim(0, 2) + 
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "red") +
  theme_minimal()
p6 <- ggplot(data.frame(value = diploid_row_means), aes(x = value)) + 
  geom_density() + 
  ggtitle("Normalized Diploid Depth") + 
  theme_minimal()
p7 <- ggplot(data.frame(value = tetraploid_row_means), aes(x = value)) + 
  geom_density() + 
  ggtitle("Normalized Tetraploid Depth") + 
  theme_minimal()
p8 <- ggplot(data.frame(value = hexaploid_row_means), aes(x = value)) + 
  geom_density() + 
  ggtitle("Normalized Hexaploid Depth") + 
  theme_minimal()
p9 <- ggplot(data.frame(value = diploid_tetraploid_ratio_cleaned), aes(x = value)) + 
  geom_density() + 
  ggtitle("Normalized Diploid - Tetraploid Ratio") + 
  xlim(0, 2) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal()
p10 <- ggplot(data.frame(value = diploid_hexaploid_ratio_cleaned), aes(x = value)) + 
  geom_density() + 
  ggtitle("Normalized Diploid - Hexaploid Ratio") + 
  xlim(0, 2) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal()
# Arrange the plots in a grid
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 2)
# Close the PDF device
dev.off()
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get a table of total depth statistics
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate statistics
total_number_chromosomes <- (16*2) + (6*4) + (4*6)
# The total number of reads reported in column three is for the whole VCF so we need to use only one of the sampled data
mean_dp_whole_data <- mean(sampled_diploid_data$total_depth)/total_number_chromosomes
median_dp_whole_data <- median(sampled_diploid_data$total_depth)/total_number_chromosomes
sd_dp_whole_data <- sd(sampled_diploid_data$total_depth)/total_number_chromosomes
min_dp_whole_data <- min(sampled_diploid_data$total_depth)/total_number_chromosomes
max_dp_whole_data <- max(sampled_diploid_data$total_depth)/total_number_chromosomes

# Create a data frame to store the results
stats_table_whole_dp <- data.frame(
  Mean = mean_dp_whole_data,
  Median = median_dp_whole_data,
  SD = sd_dp_whole_data,
  MIN = min_dp_whole_data,
  MAX = max_dp_whole_data
)
# Write the table to a CSV file
write.csv(stats_table_whole_dp, file = "results/filtered_vcf_poly/summary_statistics_depth_allVCF.csv", row.names = FALSE)
