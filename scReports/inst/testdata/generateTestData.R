library(reticulate)
library(hdf5r)
library(Matrix)

# Sample 1
# Feature types: gene expression
# Raw data dimension: 100 barcodes, 500 genes
# Filtered data dimension: 56 barcodes, 97 genes
# Metrics: gene expression
# Metric alerts: many 
#   Error: No Cells Detected
#   Warning: Low Fraction Valid Barcodes
#   Warning: Low Fraction Reads Confidently Mapped To Transcriptome
#   Warning: High Fraction of Reads Mapped Antisense to Genes
#   Warning: Low Barcode Q30 Fraction
#   Error: Low RNA Read Q30 Fraction 
#   Error: Low Sample Index Q30 Fraction 
#   Warning: Low UMI Q30 Fraction

sampleName <- "sample1"

createExample10xH5(
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
)
convertH5toH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad")
)
filterH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad"),
  outputFile = paste0("../testdata/", sampleName, "/filtered_feature_bc_matrix.h5ad")
)
values <- c(0, 54502, 1919, 66601887,"65%", "45%", "49%", "45%", "50%", 
            "72%", "95.4%", "92.4%", "4.8%", "31.1%", "56.5%", "25%", 
            "20%", "94.9%", 18391, 6628)
createExampleMetricsSummary(
  names = gexNames(),
  values = values,
  outputFile = paste0("../testdata/", sampleName, "/metrics_summary.csv")
)


# Sample 2
# Feature types: gene expression, antibody capture
# Raw data dimension: 100 barcodes, 450 genes, 50 antibodies
# Filtered data dimension: 60 barcodes, 116 genes
# Metrics: gene expression, antibody capture
# Metric alerts: none
# Other: 
#   No percent_mito column

sampleName <- "sample2"

createExample10xH5(
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  nADT = 50
)
convertH5toH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad")
)
filterH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad"),
  outputFile = paste0("../testdata/", sampleName, "/filtered_feature_bc_matrix.h5ad"),
  removePercentMito = TRUE
)
values <- c(1206, 54961, 1835, 66283828, "98.4%", "76.0%", "96.0%", "95.0%", 
            "92.9%", "95.6%","96.5%", "94.1%", "3.9%", "29.0%", "61.2%", "58.9%",
            "1.0%", "93.0%", "19608", "6233", "11372205","9429", "98.2%",
            "36.8%", "96.2%", "94.9%", "96.3%", "95.0%", "89.9%", "56.6%", "5333",
            "0.0%", "10.1%", "63.3%", "3243")
createExampleMetricsSummary(
  names = c(gexNames(), adtNames()),
  values = values,
  outputFile = paste0("../testdata/", sampleName, "/metrics_summary.csv")
)


# Sample 3
# Feature types: gene expression, crispr guide capture
# Raw data dimension: 100 barcodes, 500 genes
# Filtered data dimension: 0 barcodes, 0 features
# Metrics: gene expression, crispr guide capture
# Metric alerts: none
# Other: 
#   No cells after minimal filtering: sample should still appear in violin plots
#   Knee plot (rna) should display error message in report.

sampleName <- "sample3"

createExample10xH5(
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  nNonZero = 200,
  nBinomMu = 0.1,
  nCRISPR  = 2
)
convertH5toH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad")
)
filterH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad"),
  outputFile = paste0("../testdata/", sampleName, "/filtered_feature_bc_matrix.h5ad")
)
values <- c("11503", "81078", "5721", "932641709", "98.3%", "28.4%", "95.9%",
            "95.0%","94.0%", "95.5%", "95.2%", "92.4%", "4.7%", "17.5%", "70.1%", 
            "67.9%", "1.0%", "92.6%", "28178", "32640", "177372413", "15419",
            "96.9%", "94.0%", "96.7%", "95.3%", "95.0%", "96.7%", "72.0%", "70.7%",
            "67.4%","10395", "1.9%", "96.0%", "89.3%", "5.0%", "486")
createExampleMetricsSummary(
  names = c(gexNames(), crisprNames()),
  values = values,
  outputFile = paste0("../testdata/", sampleName, "/metrics_summary.csv")
)


# Sample 4
# Feature types: gene expression, antibody capture
# Raw data dimension: 100 barcodes, 400 genes, 100 antibodies
# Filtered data dimension: 53 barcodes, 123 genes
# Metric alerts: many
#   Warning: Low Number of Cells Detected
#   Error: Low Fraction Valid Barcodes
#   Error: Low Fraction Reads Confidently Mapped To Transcriptome
#   Error: High Fraction of Reads Mapped Antisense to Genes
#   Error: Low Barcode Q30 Fraction 
#   Warning: Low RNA Read Q30 Fraction 
#   Warning: Low Sample Index Q30 Fraction 
#   Error: Low UMI Q30 Fraction 
# Other: 
#   Knee plot (adt) should display error message in report.

sampleName <- "sample4"

createExample10xH5(
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  nADT = 100
)
convertH5toH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad")
)
filterH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad"),
  outputFile = paste0("../testdata/", sampleName, "/filtered_feature_bc_matrix.h5ad")
)
values <- c(90, 63617, 1930, 45359427, "40%", "30%", "20%", "59%", "76%",
            "45%", "94.9%", "15%", "5.5%", "31.8%", "54.4%", "15%", "40%",
            "87.6%", 17467, 6416, 12606650, 17681, "98.5%", "48.2%", "96.6%",
            "94.4%", "90.8%", "96.0%", "90.2%", "56.4%", 9974, "0.0%", 
            "9.8%", "63.1%", 4404)
createExampleMetricsSummary(
  names = c(gexNames(), adtNames()),
  values = values,
  outputFile = paste0("../testdata/", sampleName, "/metrics_summary.csv")
)


# Sample 5
# Feature types: gene expression (cell ranger 2)
# Raw data dimension: 100 barcodes, 500 genes
# Filtered data dimension:  56 barcodes, 97 genes
# Metric alerts: none

sampleName <- "sample5"

createExample10xH5(
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  groupName = "mm10",
  version = 2
)
convertH5toH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5"),
  outputFile = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad")
)
filterH5AD(
  filePath = paste0("../testdata/", sampleName, "/raw_feature_bc_matrix.h5ad"),
  outputFile = paste0("../testdata/", sampleName, "/filtered_feature_bc_matrix.h5ad")
)
values <- c(931, "56718", "2798", "52805264", "98.1%", "57.3%", "97.3%",
            "86.0%", "97.3%", "97.7%", "94.1%", "89.2%", "3.3%", "20.7%",
            "65.1%", "62.3%", "1.3%", "78.3%", "16152", "8545")
createExampleMetricsSummary(
  names = gexNames(),
  values = values,
  outputFile = paste0("../testdata/", sampleName, "/metrics_summary.csv")
)

# Meta

createMetaFile(
  nSamples = 5,
  outputFile = "../testdata/meta/meta.tsv"
)


# Sample 6
# Feature types: gene expression
# Raw data dimension: 0 barcodes, 0 genes
# Filtered data dimension:  0 barcodes, 0 genes
# Metric alerts: none

# Does 10x return and h5 file when no cells are present?



