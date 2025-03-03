#!/usr/bin/env Rscript

# run_pipeline.R - Main script to run the entire viral protein database pipeline
# This script orchestrates the complete pipeline from downloading polyproteins
# to BLAST annotation of mature proteins

# Load all pipeline components
source("01_config.R")
source("02_download_polyproteins.R")
source("03_extract_mature_proteins.R")
source("04_blast_annotation.R")

#' Run the complete viral protein database pipeline
#' 
#' @param max_records Maximum number of records to process (NULL for all)
#' @param email User email for NCBI (NULL to prompt)
#' @param api_key NCBI API key (NULL for none or to prompt)
#' @return List of output files generated
#' @export
run_viral_protein_pipeline <- function(max_records = NULL, email = NULL, api_key = NULL) {
  cat("\n=== Viral Protein Database Pipeline ===\n\n")
  
  # 1. Set up credentials
  credentials <- NULL
  if(!is.null(email)) {
    credentials <- list(email = email, api_key = api_key)
    saveRDS(credentials, "viral_db_credentials.rds")
    cat("Using provided credentials for NCBI access.\n")
  } else if(file.exists("viral_db_credentials.rds")) {
    credentials <- readRDS("viral_db_credentials.rds")
    cat("Using stored credentials for NCBI access.\n")
  } else {
    cat("No credentials provided. Requesting from user...\n")
    credentials <- get_ncbi_credentials()
  }
  
  # Create date-stamped filenames
  today <- format(Sys.Date(), "%Y_%m_%d")
  polyprotein_file <- paste0(today, "_viral_polyprotein.fasta")
  mature_protein_file <- paste0(today, "_viral_polyprotein_mature_proteins.fasta")
  missing_annotation_file <- paste0(today, "_viral_polyprotein_missing_annotation.fasta")
  blast_result_file <- paste0(today, "_blastx_results.txt")
  
  # 2. Download polyproteins
  cat("\n=== Step 1: Downloading Viral Polyproteins ===\n\n")
  tryCatch({
    polyprotein_file <- download_viral_polyproteins(
      credentials = credentials,
      max_records = max_records,
      output_file = polyprotein_file
    )
    cat("Polyprotein download completed successfully.\n")
  }, error = function(e) {
    cat("ERROR in polyprotein download:", conditionMessage(e), "\n")
    stop("Pipeline failed at step 1.")
  })
  
  # 3. Extract mature proteins
  cat("\n=== Step 2: Extracting Mature Proteins ===\n\n")
  tryCatch({
    output_files <- extract_mature_proteins(
      credentials = credentials,
      max_records = max_records,
      output_file = mature_protein_file,
      missing_file = missing_annotation_file
    )
    cat("Mature protein extraction completed successfully.\n")
  }, error = function(e) {
    cat("ERROR in mature protein extraction:", conditionMessage(e), "\n")
    stop("Pipeline failed at step 2.")
  })
  
  # 4. Run BLAST annotation
  cat("\n=== Step 3: BLAST Annotation of Unannotated Polyproteins ===\n\n")
  blast_success <- FALSE
  tryCatch({
    blast_success <- run_blast_annotation(
      mature_proteins_file = mature_protein_file,
      missing_annotations_file = missing_annotation_file,
      taxonomy_map_file = "viral_taxonomy_map.rds",
      output_blast = blast_result_file
    )
    if(blast_success) {
      cat("BLAST annotation completed successfully.\n")
    } else {
      cat("BLAST annotation completed with warnings.\n")
    }
  }, error = function(e) {
    cat("ERROR in BLAST annotation:", conditionMessage(e), "\n")
    cat("Pipeline will continue despite BLAST errors.\n")
  })
  
  # 5. Wrap up and return results
  cat("\n=== Pipeline Complete ===\n\n")
  cat("The following files were generated:\n")
  cat("1. Viral polyproteins:", polyprotein_file, "\n")
  cat("2. Mature proteins:", mature_protein_file, "\n")
  cat("3. Unannotated polyproteins:", missing_annotation_file, "\n")
  if(blast_success) {
    cat("4. BLAST results:", blast_result_file, "\n")
  }
  
  return(list(
    polyproteins = polyprotein_file,
    mature_proteins = mature_protein_file,
    missing_annotations = missing_annotation_file,
    blast_results = if(blast_success) blast_result_file else NULL
  ))
}

# If this script is run directly, run the complete pipeline
if(!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  max_records <- NULL
  email <- NULL
  api_key <- NULL
  
  # Simple argument parsing
  if(length(args) >= 1) {
    if(!is.na(as.numeric(args[1]))) {
      max_records <- as.numeric(args[1])
    } else {
      email <- args[1]
    }
  }
  
  if(length(args) >= 2) {
    if(is.null(email)) {
      email <- args[1]
    }
    api_key <- args[2]
  }
  
  # Run the pipeline
  cat("Starting viral protein database pipeline...\n")
  if(!is.null(max_records)) {
    cat("Limiting to", max_records, "records for testing.\n")
  }
  
  results <- run_viral_protein_pipeline(max_records, email, api_key)
  
  cat("\nPipeline execution complete.\n")
  cat("All output files are ready for use.\n")
}