#' ---
#'
#' This script extracts the GeneScoreMatrix from an ArchR project for each autosome (chr1–chr22),
#' and saves the matrix, genes, and barcodes in Matrix Market format for downstream analysis.
#' 
#' Usage: Run in an R environment with access to the appropriate ArchR project.
#'
#' Output: Creates a directory per chromosome with:
#'         - matrix.mtx (gene scores)
#'         - genes.tsv (gene names)
#'         - barcodes.tsv (cell names)
#' 
#' Author: Hessel Dijkstra
#' ---

library(ArchR)   
library(Matrix)  

# Set working directory to the location where ArchR projects are stored
setwd("../../../archr_preprocess_samples/")

# Configure ArchR
addArchRThreads(threads = 4)
addArchRGenome("hg38")
h5disableFileLocking()

# Load the preprocessed ArchR project
proj <- loadArchRProject("Project_common_atac_rna_cells")

# Define output directory for the exported data
out_dir <- "../../../sc-scoring-student-project"

# Specify the chromosomes to process (autosomes only)
chromosomes <- paste0("chr", 1:22)

# Loop over each chromosome and extract gene scores
for (chr in chromosomes) {
  message("Processing chromosome: ", chr)
  
  # Extract the GeneScoreMatrix for the current chromosome
  gs <- getMatrixFromProject(
    proj,
    useMatrix = "GeneScoreMatrix",  # Type of matrix to extract
    useSeqnames = chr,              # Chromosome to filter by
    logFile = file.path(out_dir, "archr_log.txt")  # Log file for ArchR
  )
  
  # Access the sparse gene score matrix
  mat <- assays(gs)$GeneScoreMatrix
  
  # Get gene names and cell barcodes
  genes <- rowData(gs)$name
  cells <- colnames(mat)
  
  # Create chromosome-specific output directory
  chr_dir <- file.path(out_dir, sprintf("gene_score_batch_%s", chr))
  dir.create(chr_dir, showWarnings = FALSE)
  
  # Write matrix, genes, and barcodes to files in Matrix Market format
  writeMM(mat, file.path(chr_dir, "matrix.mtx"))  # Sparse matrix file
  write.table(genes, file.path(chr_dir, "genes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)  # Gene names
  write.table(cells, file.path(chr_dir, "barcodes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)  # Cell barcodes
  
  print("files written")
}
