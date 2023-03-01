sampleNames <- c(
  "sample1",
  "sample2",
  "sample3",
  "sample4",
  "sample5"
)
sampleCRQC <- c(
  "../testdata/sample1/metrics_summary.csv", 
  "../testdata/sample2/metrics_summary.csv", 
  "../testdata/sample3/metrics_summary.csv",
  "../testdata/sample4/metrics_summary.csv",
  "../testdata/sample5/metrics_summary.csv"
)
metricsReference <- data.frame(
  "SampleID" = c("sample1", "sample2", "sample3", "sample4", "sample5"),
  "Estimated Number of Cells" = c(0, 1206, 11503, 90, 931),
  "Mean Reads per Cell" = c(54502, 54961, 81078, 63617, 56718),
  "Median Genes per Cell" = c(1919, 1835, 5721, 1930, 2798),
  "Number of Reads" = c(66601887, 66283828, 932641709, 45359427, 52805264),
  "Valid Barcodes" = c(65, 98.4, 98.3, 40.0, 98.1),
  "Sequencing Saturation" = c(45, 76, 28.4, 30.0, 57.3),
  "Q30 Bases in Barcode" = c(49, 96, 95.9, 20.0, 97.3),
  "Q30 Bases in RNA Read" = c(45, 95, 95, 59, 86),
  "Q30 Bases in Sample Index" = c(50, 92.9, 94, 76, 97.3),
  "Q30 Bases in UMI" = c(72, 95.6, 95.5, 45, 97.7),
  "Reads Mapped to Genome" = c(95.4, 96.5, 95.2, 94.9, 94.1),
  "Reads Mapped Confidently to Genome" = c(92.4, 94.1, 92.4, 15, 89.2),
  "Reads Mapped Confidently to Intergenic Regions" = c(4.8, 3.9, 4.7, 5.5, 3.3),
  "Reads Mapped Confidently to Intronic Regions" = c(31.1, 29, 17.5, 31.8, 20.7),
  "Reads Mapped Confidently to Exonic Regions" = c(56.5, 61.2, 70.1, 54.4, 65.1),
  "Reads Mapped Confidently to Transcriptome" = c(25, 58.9, 67.9, 15, 62.3),
  "Reads Mapped Antisense to Gene" = c(20, 1, 1, 40, 1.3),
  "Fraction Reads in Cells" = c(94.9, 93, 92.6, 87.6, 78.3),
  "Total Genes Detected" = c(18391, 19608, 28178, 17467, 16152),                                                         
  "Median UMI Counts per Cell" = c(6628, 6233, 32640, 6416, 8545),
  "Antibody: Number of Reads" = c(NA, 11372205, NA, 12606650, NA),                                                    
  "Antibody: Mean Reads per Cell" = c(NA, 9429, NA, 17681, NA),
  "Antibody: Valid Barcodes" = c(NA, 98.2, NA, 98.5, NA),                                                      
  "Antibody: Sequencing Saturation" = c(NA, 36.8, NA, 48.2, NA),
  "Antibody: Q30 Bases in Barcode" = c(NA, 96.2, NA, 96.6, NA),                                               
  "Antibody: Q30 Bases in Antibody Read" = c(NA, 94.9, NA, 94.4, NA),
  "Antibody: Q30 Bases in UMI" = c(NA, 96.3, NA, 90.8, NA), 
  "Antibody: Q30 Bases in Sample Index" = c(NA, 95, NA, 96, NA),
  "Antibody: Fraction Antibody Reads" = c(NA, 89.9, NA, 90.2, NA),
  "Antibody: Fraction Antibody Reads Usable"  = c(NA, 56.6, NA, 56.4, NA),                                    
  "Antibody: Antibody Reads Usable per Cell" = c(NA, 5333, NA, 9974, NA),
  "Antibody: Fraction Reads in Barcodes with High UMI Counts" = c(NA, 0, NA, 0, NA),                   
  "Antibody: Fraction Unrecognized Antibody" = c(NA, 10.1, NA, 9.8, NA),
  "Antibody: Antibody Reads in Cells" = c(NA, 63.3, NA, 63.1, NA),                                          
  "Antibody: Median UMIs per Cell (summed over all recognized antibody barcodes)" = c(NA, 3243, NA, 4404, NA),
  "CRISPR: Number of Reads" = c(NA, NA, 177372413, NA, NA),
  "CRISPR: Mean Reads per Cell" = c(NA, NA, 15419, NA, NA),                                                
  "CRISPR: Valid Barcodes" = c(NA, NA, 96.9, NA, NA),
  "CRISPR: Sequencing Saturation" = c(NA, NA, 94, NA, NA),                                               
  "CRISPR: Q30 Bases in Barcode" = c(NA, NA, 96.7, NA, NA),
  "CRISPR: Q30 Bases in CRISPR Read" = c(NA, NA, 95.3, NA, NA),
  "CRISPR: Q30 Bases in Sample Index" =  c(NA, NA, 95, NA, NA),
  "CRISPR: Q30 Bases in UMI" = c(NA, NA, 96.7, NA, NA),
  "CRISPR: Fraction Reads with Putative Protospacer Sequence" = c(NA, NA, 72, NA, NA),                  
  "CRISPR: Fraction Guide Reads" = c(NA, NA, 70.7, NA, NA),
  "CRISPR: Fraction Guide Reads Usable"  = c(NA, NA, 67.4, NA, NA),                                         
  "CRISPR: Guide Reads Usable per Cell" = c(NA, NA, 10395, NA, NA),
  "CRISPR: Fraction Protospacer Not Recognized" = c(NA, NA, 1.9, NA, NA),                                
  "CRISPR: Guide Reads in Cells" = c(NA, NA, 96, NA, NA),
  "CRISPR: Cells with 1 or more protospacers detected" = c(NA, NA, 89.3, NA, NA),                          
  "CRISPR: Cells with 2 or more protospacers detected" = c(NA, NA, 5, NA, NA),
  "CRISPR: Median UMIs per Cell (summed over all recognized protospacers)" = c(NA, NA, 486, NA, NA)
)


metrics <- readCRQC(sampleNames = sampleNames, 
                    sampleCRQC = sampleCRQC)

testthat::test_that(
  "a data.frame is returned that is equal to the reference data.frame", {
    
  testthat::expect_s3_class(object = metrics, class = "data.frame")
    
  testthat::expect_equal(object = metrics, 
                         expected = metricsReference, 
                         ignore_attr = TRUE)
  
  }
)

