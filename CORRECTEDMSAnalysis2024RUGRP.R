#Package Installation
install.packages("BiocManager") #BioConductor Manager package
install.packages(c("MSnbase","limma","DESeq2","ggplot2","dplyr","UpSetR","ComplexHeatmap"))
install.packages("tidyverse")
##Library
library(BiocManager)
library(MSnbase)
library(limma)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(FragPipeAnalystR)
library(tidyverse)


#Import Data
## Data important for proteomics (and for this specific research project) would have to include the protein group matrix files from DIA outputs and combined_protein files from the DDA.
## Here we use tidyverse to read the tsv files into a table that can be read easily and we are able to extract columns or specific cells from said table.
### DDA outputs - Exploris
EXPLORIS_DDA_combined_ion <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DDA_analysis_output/combined_ion.tsv") 
EXPLORIS_DDA_combined_peptide <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DDA_analysis_output/combined_peptide.tsv")
EXPLORIS_DDA_combined_protein <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DDA_analysis_output/combined_protein.tsv")
EXPLORIS_DDA_combined_modified_peptide <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DDA_analysis_output/combined_modified_peptide.tsv")

### DDA outputs - TimsTOF (15cm)
TIMSTOF15CM_DDA_combined_ion <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DDA_15cm_analysis_output/combined_ion.tsv")
TIMSTOF15CM_DDA_combined_peptide <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DDA_15cm_analysis_output/combined_peptide.tsv")
TIMSTOF15CM_DDA_combined_protein <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DDA_15cm_analysis_output/combined_protein.tsv")
TIMSTOF15CM_DDA_combined_modified_peptide <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DDA_15cm_analysis_output/combined_modified_peptide.tsv")

### DIA outputs - Exploris 
EXPLORIS_DIA_STAT_REPORT <-read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DIA_analysis_output/diann-output/report.stats.tsv")
EXPLORIS_DIA_PROTEIN_GROUP <-read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DIA_analysis_output/diann-output/report.pg_matrix.tsv")

### DIA outputs - TimsTOF (15cm)
TIMSTOF15CM_DIA_STAT_REPORT <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_15cm_analysis_output/diann-output/report.stats.tsv")
TIMSTOF15CM_DIA_PROTEIN_GROUP <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_15cm_analysis_output/diann-output/report.pg_matrix.tsv")
TIMSTOF5CM_DIA_EXPERIMENT_ANNOTATION <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_5cm_analysis_output/experiment_annotation.tsv")

### DIA outputs - TImsTOF (5cm)
TIMSTOF5CM_DIA_STAT_REPORT <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_5cm_analysis_output/diann-output/report.stats.tsv")
TIMSTOF5CM_DIA_PROTEIN_GROUP <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_5cm_analysis_output/diann-output/report.pg_matrix.tsv")
TIMSTOF15CM_DIA_EXPERIMENT_ANNOTATION <- read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_15cm_analysis_output/experiment_annotation.tsv")

#Quality Control (QC)
## Histograms for the precursors identified
ggplot(EXPLORIS_DIA_STAT_REPORT, aes(x = Precursors.Identified)) +
  geom_histogram(binwidth = 50, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Precursors Identified",
       x = "Number of Precursors Identified",
       y = "Frequency") +
  theme_minimal()

ggplot(TIMSTOF5CM_DIA_STAT_REPORT, aes(x = Precursors.Identified)) +
  geom_histogram(binwidth = 50, fill = "yellow", color = "black") +
  labs(title = "Histogram of Precursors Identified",
       x = "Number of Precursors Identified",
       y = "Frequency") +
  theme_minimal()

ggplot(TIMSTOF15CM_DIA_STAT_REPORT, aes(x = Precursors.Identified)) +
  geom_histogram(binwidth = 50, fill = "red", color = "black") +
  labs(title = "Histogram of Precursors Identified",
       x = "Number of Precursors Identified",
       y = "Frequency") +
  theme_minimal()

### We can also compare the three 
precursors_exploris <- EXPLORIS_DIA_STAT_REPORT$Precursors.Identified
precursors_timstof5cm <-TIMSTOF5CM_DIA_STAT_REPORT$Precursors.Identified
precursors_timstof15cm <-TIMSTOF15CM_DIA_STAT_REPORT$Precursors.Identified

#### calculating summary statistics
summary_stats <- data.frame(
  Platform = c("Exploris", "TimsTOF 15cm", "TimsTOF 5cm"),
  Mean = c(mean(precursors_exploris), mean(precursors_timstof15cm), mean(precursors_timstof5cm)),
  Median = c(median(precursors_exploris), median(precursors_timstof15cm), median(precursors_timstof5cm)),
  SD = c(sd(precursors_exploris), sd(precursors_timstof15cm), sd(precursors_timstof5cm))
)

ggplot(summary_stats, aes(x = Platform, y = Mean)) +
  geom_bar(stat = "identity", fill = "lightblue", width = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "black") +
  geom_point(aes(y = Median), color = "red", size = 3) +
  labs(title = "Summary Statistics for Precursors Identified",
       x = "Platform",
       y = "Mean",
       caption = "Error bars represent standard deviation") +
  theme_minimal()

# Analysis
##Column Comparison between DIA Timstof 15cm and 5cm
### Doing T test on number of proteins identified
protein_number_5cm <- TIMSTOF5CM_DIA_STAT_REPORT$Proteins.Identified
protein_number_15cm <- TIMSTOF15CM_DIA_STAT_REPORT$Proteins.Identified
column_t_test <-t.test(protein_number_5cm, protein_number_15cm )

se15cm <- make_se_from_files("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_15cm_analysis_output/diann-output/report.pg_matrix.tsv",
                             "C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_15cm_analysis_output/experiment_annotation.tsv", level="protein", type = "DIA" )


se5cm <- make_se_from_files("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_5cm_analysis_output/diann-output/report.pg_matrix.tsv",
                            "C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/TimsTOF_DIA_5cm_analysis_output/experiment_annotation.tsv", level="protein", type = "DIA" )
plot_pca(se15cm)
plot_pca(se5cm)

plot_correlation_heatmap(se15cm)
plot_correlation_heatmap(se5cm)

{boxplot(protein_number_5cm, protein_number_15cm, names = c("5cm column", "15 cm column"),
         xlab = "Dataset", ylab = "Number of Proteins Identified",
         col = c("skyblue", "lightgreen"))
  if (column_t_test$p.value < 0.05) {
    text(1.5, max(max(protein_number_5cm), max(protein_number_15cm)), "*", cex = 2)
  }
  
  
  p_value_text <- paste("p-value:", round(column_t_test$p.value, 30))
  alpha_text <- paste("alpha:", 0.05)
  significance_text <- ifelse(column_t_test$p.value < 0.05, "Significant", "Not Significant")
  
  text(1, max(max(protein_number_5cm), max(protein_number_15cm)) - 100, p_value_text, adj = 0, cex = 0.8)
  text(1, max(max(protein_number_5cm), max(protein_number_15cm)) - 200, alpha_text, adj = 0, cex = 0.8)
  text(1, max(max(protein_number_5cm), max(protein_number_15cm)) - 300, significance_text, adj = 0, cex = 0.8)}

#### Exploratory Analysis in time/sample amount/protein count
TIMSTOF5CM_DIA_STAT_REPORT$Sample_amount_ng <- 10
TIMSTOF15CM_DIA_STAT_REPORT$Sample_amount_ng <-10
TIMSTOF5CM_DIA_STAT_REPORT$Run_time_minutes <- 15
TIMSTOF15CM_DIA_STAT_REPORT$Run_time_minutes <- 32

barplot(c(sum(TIMSTOF5CM_DIA_STAT_REPORT$Proteins.Identified), sum(TIMSTOF15CM_DIA_STAT_REPORT$Proteins.Identified)),
        names.arg = c("5cm Column", "15cm Column"),
        xlab = "Column Length",
        ylab = "Total Proteins Identified",
        main = "Comparison of Total Proteins Identified")


#### Efficiency Metrics  
TIMSTOF5CM_DIA_STAT_REPORT$Proteins_per_minute <- TIMSTOF5CM_DIA_STAT_REPORT$Proteins.Identified / 15  # Assuming 15 minutes for 5cm column
TIMSTOF5CM_DIA_STAT_REPORT$Proteins_per_ng <- TIMSTOF5CM_DIA_STAT_REPORT$Proteins.Identified / 50     # Assuming 50ng of sample for 5cm column
TIMSTOF15CM_DIA_STAT_REPORT$Proteins_per_minute <- TIMSTOF15CM_DIA_STAT_REPORT$Proteins.Identified / 32  # Assuming 32 minutes for 15cm column
TIMSTOF15CM_DIA_STAT_REPORT$Proteins_per_ng <- TIMSTOF15CM_DIA_STAT_REPORT$Proteins.Identified / 50     # Assuming 50ng of sample for 15cm column

{plot(TIMSTOF5CM_DIA_STAT_REPORT$Proteins_per_minute, TIMSTOF15CM_DIA_STAT_REPORT$Proteins_per_minute,
     xlab = "Proteins per minute (5cm)", ylab = "Proteins per minute (15cm)",
     main = "Proteins per minute Comparison")
abline(0, 1, col = "red")}


####PRINTING NUMBER OF PROTEINS IN EACH MACHINE
num_proteins_EXPLORIS_DDA <- nrow(EXPLORIS_DDA_combined_protein)
num_proteins_TIMSTOF15CM_DDA <-nrow(TIMSTOF15CM_DDA_combined_protein)
cat("Number of identified proteins in Exploris DDA:", num_proteins_EXPLORIS_DDA, "\n")
cat("Number of identified proteins in TimsTOF DDA:", num_proteins_TIMSTOF15CM_DDA, "\n")

TIMSTOF15CM_DDA_abundances <- TIMSTOF15CM_DDA_combined_peptide[, grep("^[TV]([1-3]_[1-3])? Intensity$", colnames(TIMSTOF15CM_DDA_combined_peptide))]
EXPLORIS_DDA_abundances <- EXPLORIS_DDA_combined_peptide[, grep("^[TV](_[1-3])? Intensity$", colnames(EXPLORIS_DDA_combined_peptide))]

print(colnames(EXPLORIS_DDA_combined_peptide))
print(colnames(TIMSTOF15CM_DDA_combined_peptide))
print(TIMSTOF15CM_DIA_PROTEIN_GROUP)

#Trans-platform analysis between DDA exploris vs timstof
EXPLORIS_DIA_STAT_REPORT <-read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DIA_analysis_output/diann-output/report.stats.tsv")
EXPLORIS_DIA_PROTEIN_GROUP <-read_tsv("C:/Users/Umar/Desktop/R_Studio/RP_RUG_2024/Fragpipe_output/Exploris_DIA_analysis_output/diann-output/report.pg_matrix.tsv")

#Trans-platform analysis between DDA exploris vs timstof
EXPLORIS_DDA_TOTAL_PROTEIN_T <- EXPLORIS_DDA_combined_protein$`Protein ID`  #extracting protein ID from exploris
TIMSTOF15CM_DDA_TOTAL_PROTEIN_T <- TIMSTOF15CM_DDA_combined_protein$`Protein ID` #extracting protein ID from timstof 15cm
EXPLORIS_DIA_TOTAL_PROTEIN_T <- EXPLORIS_DIA_PROTEIN_GROUP$Protein.Ids #extract protein id
TIMSTOF15CM_DIA_TOTAL_PROTEIN_T <- TIMSTOF15CM_DIA_PROTEIN_GROUP$Protein.Ids
TIMSTOF5CM_DIA_TOTAL_PROTEIN_T <- TIMSTOF5CM_DIA_PROTEIN_GROUP$Protein.Ids


UPSET_INPUT <- list(Exploris_DDA=EXPLORIS_DDA_TOTAL_PROTEIN_T, TimsTOF15cm_DDA=TIMSTOF15CM_DDA_TOTAL_PROTEIN_T, TimsTOF15cm_DIA=TIMSTOF15CM_DIA_TOTAL_PROTEIN_T,
                    TimsTOF5cm_DIA=TIMSTOF5CM_DIA_TOTAL_PROTEIN_T, Exploris_DIA=EXPLORIS_DIA_TOTAL_PROTEIN_T)

list_to_matrix(UPSET_INPUT)

upset(fromList(UPSET_INPUT), order.by = "freq")

{ggplot(summary_stats, aes(x = Platform, y = Mean)) +
  geom_bar(stat = "identity", fill = "lightblue", width = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "black") +
  geom_point(aes(y = Median), color = "red", size = 3) +
  labs(title = "Summary Statistics for Precursors Identified",
       x = "Platform",
       y = "Mean",
       caption = "Error bars represent standard deviation") +
  theme_minimal()}


#MAKE SURE TO ADD DISTINCT() TO NORMALIZE DATA TO REMOVE DUPLICATES.
