#!/usr/bin/env Rscript

# 01_config.R - Configuration settings for viral protein database pipeline
# This script sets up the environment and user credentials for the pipeline

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Required packages
required_packages <- c("rentrez", "Biostrings", "stringr", "seqinr", "dplyr")

# Install missing packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    if(pkg == "Biostrings") {
      if (!require("BiocManager")) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Get NCBI credentials from user
#' 
#' @return A list containing email and API key
#' @export
get_ncbi_credentials <- function() {
  cat("NCBI E-utilities requires your email for better service.\n")
  
  # Get email with validation
  email <- ""
  while(email == "") {
    email <- readline(prompt = "Enter your email address: ")
    if(email == "" || !grepl("^[[:alnum:]_.+-]+@[[:alnum:]_.-]+\\.[[:alpha:]]+$", email)) {
      cat("Please enter a valid email address.\n")
      email <- ""
    }
  }
  
  # Get optional API key
  use_api_key <- readline(prompt = "Do you have an NCBI API key? (y/n): ")
  api_key <- NULL
  
  if(tolower(substr(use_api_key, 1, 1)) == "y") {
    api_key <- readline(prompt = "Enter your NCBI API key: ")
  }
  
  # Create credentials list
  credentials <- list(
    email = email,
    api_key = api_key
  )
  
  # Save credentials to temporary file for use across scripts
  saveRDS(credentials, "viral_db_credentials.rds")
  
  cat("Credentials saved successfully.\n")
  return(credentials)
}

#' Create standardized date-stamped filename
#' 
#' @param base_name Base name for the file
#' @param extension File extension (default: "fasta")
#' @return Filename with date stamp
#' @export
get_dated_filename <- function(base_name, extension = "fasta") {
  today <- format(Sys.Date(), "%Y_%m_%d")
  return(paste0(today, "_", base_name, ".", extension))
}

#' Format taxonomy with placeholders
#' 
#' @param taxonomy Taxonomy string from NCBI
#' @return Vector of standardized taxonomy ranks with placeholders
#' @export
format_taxonomy <- function(taxonomy) {
  ranks <- c("phylum", "class", "order", "family", "genus", "species")
  result <- rep("Unknown", length(ranks))
  names(result) <- ranks
  
  if(is.null(taxonomy) || taxonomy == "") {
    return(result)
  }
  
  # Split taxonomy by semicolons and clean
  tax_parts <- strsplit(taxonomy, ";")[[1]]
  tax_parts <- trimws(tax_parts)
  
  # Try to match each rank
  for(i in seq_along(ranks)) {
    rank <- ranks[i]
    matched <- FALSE
    
    # Direct pattern matching
    pattern <- paste0(rank, ":")
    matches <- grepl(pattern, tax_parts, ignore.case = TRUE)
    
    if(any(matches)) {
      idx <- which(matches)[1]
      result[i] <- gsub(".*:", "", tax_parts[idx])
      matched <- TRUE
    }
    
    # Special matching for common viral taxonomy patterns
    if(!matched) {
      if(rank == "family" && any(grepl("viridae$", tax_parts, ignore.case = TRUE))) {
        idx <- grep("viridae$", tax_parts, ignore.case = TRUE)[1]
        result[i] <- tax_parts[idx]
        matched <- TRUE
      } else if(rank == "genus" && any(grepl("virus$", tax_parts, ignore.case = TRUE)) && 
                !any(grepl("viridae$|virales$", tax_parts, ignore.case = TRUE))) {
        idx <- grep("virus$", tax_parts, ignore.case = TRUE)[1]
        result[i] <- tax_parts[idx]
        matched <- TRUE
      }
    }
    
    # Create placeholder based on higher taxonomy if not matched
    if(!matched && i > 1 && result[i-1] != "Unknown") {
      result[i] <- paste0(result[i-1], "_no_", rank)
    }
  }
  
  # Replace spaces with underscores
  result <- gsub(" ", "_", result)
  
  return(result)
}

#' Create standardized FASTA header
#' 
#' @param accession Accession number
#' @param taxonomy Taxonomy string
#' @param organism Organism name
#' @param protein_name Protein name
#' @return Formatted FASTA header
#' @export
create_fasta_header <- function(accession, taxonomy, organism, protein_name) {
  # Format taxonomy
  tax_parts <- format_taxonomy(taxonomy)
  
  # Clean organism and protein names (replace spaces with underscores)
  organism <- gsub(" ", "_", organism)
  protein_name <- gsub(" ", "_", protein_name)
  
  # Create header
  header <- paste0(
    ">", accession, "|Viruses|",
    paste(tax_parts, collapse = "|"), "|",
    organism, "|", protein_name
  )
  
  return(header)
}

# If this script is run directly (not sourced), ask for credentials
if(!interactive()) {
  cat("Initializing viral protein database pipeline...\n")
  credentials <- get_ncbi_credentials()
  cat("Configuration complete. Ready to proceed with pipeline.\n")
}