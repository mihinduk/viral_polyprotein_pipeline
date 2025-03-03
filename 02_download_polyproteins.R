#!/usr/bin/env Rscript

# 02_download_polyproteins.R - Download viral polyproteins from NCBI
# This script downloads complete viral polyproteins and saves them to a FASTA file

# Load configuration and helper functions
source("01_config.R")

#' Download all complete viral polyproteins
#' 
#' @param credentials NCBI credentials (email and API key)
#' @param max_records Maximum number of records to download (NULL for all)
#' @param output_file Output FASTA filename (NULL for auto-generated)
#' @return Path to output file with downloaded polyproteins
#' @export
download_viral_polyproteins <- function(credentials = NULL, max_records = NULL, output_file = NULL) {
  # Load credentials if not provided
  if(is.null(credentials)) {
    if(file.exists("viral_db_credentials.rds")) {
      credentials <- readRDS("viral_db_credentials.rds")
    } else {
      credentials <- get_ncbi_credentials()
    }
  }
  
  # Set up rentrez with credentials
  if(!is.null(credentials$api_key) && credentials$api_key != "") {
    set_entrez_key(credentials$api_key)
    cat("Using API key for faster download rates.\n")
  }
  
  # Set output filename
  if(is.null(output_file)) {
    output_file <- get_dated_filename("viral_polyprotein")
  }
  
  # Create or clear output file
  file.create(output_file)
  cat("Output will be saved to:", output_file, "\n")
  
  # Search term for viral polyproteins
  search_term <- "(virus[organism] OR viruses[organism]) AND polyprotein[protein] AND complete[title]"
  
  # Perform initial search to get count
  cat("Searching NCBI for viral polyproteins...\n")
  search_results <- entrez_search(db = "nuccore", term = search_term, retmax = 0)
  total_count <- search_results$count
  
  if(total_count == 0) {
    cat("No records found matching the search criteria.\n")
    return(output_file)
  }
  
  cat("Found", total_count, "viral polyprotein records.\n")
  
  # Limit records if requested
  if(!is.null(max_records) && max_records < total_count) {
    cat("Limiting download to", max_records, "records.\n")
    total_count <- max_records
  }
  
  # Initialize counters
  polyprotein_count <- 0
  batch_size <- 100
  skipped_count <- 0
  
  # Process in batches
    for(i in seq(1, total_count, by = batch_size)) {
      batch_end <- min(i + batch_size - 1, total_count)
      cat(paste0("Processing records ", i, " to ", batch_end, " of ", total_count, "...\n"))

      # Get batch of IDs
      this_batch_size <- min(batch_size, total_count - (i - 1))
      cat("Requesting", this_batch_size, "records in this batch...\n")
      search_batch <- entrez_search(
        db = "nuccore",
        term = search_term,
        retmax = this_batch_size,
        retstart = i - 1,
        use_history = TRUE
      )
      cat("Got", length(search_batch$ids), "record IDs\n")
    
    # Fetch records
    records <- entrez_fetch(
      db = "nuccore", 
      rettype = "gb", 
      retmode = "text",
      web_history = search_batch$web_history,
      parsed = FALSE
    )
    
    # Split into individual GenBank records
    gb_records <- strsplit(records, "LOCUS")[[1]][-1]
    
    if(length(gb_records) == 0) {
      cat("No records found in this batch.\n")
      next
    }
    
    cat("Processing", length(gb_records), "records in this batch.\n")

    # Add debugging
  cat("First 100 characters of records data:\n")
  cat(substr(records, 1, 100), "\n")
  cat("Number of records after splitting:", length(gb_records), "\n")
    
    # Process each record
    for(record_text in gb_records) {
      # Add LOCUS back for proper parsing
      record_text <- paste0("LOCUS", record_text)
      
      # Extract basic info
      accession <- gsub(".*ACCESSION\\s+([^\\s]+).*", "\\1", record_text, perl = TRUE)
      if(accession == record_text) accession <- "unknown"
      
      organism_match <- regexpr("ORGANISM\\s+([^\\n]+)", record_text, perl = TRUE)
      organism <- if(organism_match > 0) {
        trimws(gsub("ORGANISM\\s+", "", regmatches(record_text, organism_match)))
      } else {
        "unknown"
      }
      
      # Extract taxonomy
      taxonomy_block <- gsub(".*ORGANISM[^\\n]*\\n\\s*([^/]*?)\\n/.*", "\\1", record_text, perl = TRUE)
      taxonomy <- gsub("\\s+", " ", taxonomy_block)
      
      # Find CDS features that contain "polyprotein"
      cds_pattern <- "/CDS\\s+([^/]+).*?/product=\"([^\"]*polyprotein[^\"]*)\".*?(/protein_id=\"([^\"]+)\")?.*?\\n"
      cds_matches <- gregexpr(cds_pattern, record_text, perl = TRUE, ignore.case = TRUE)
      
      # Process CDS features
      if(cds_matches[[1]][1] != -1) {
        cds_texts <- regmatches(record_text, cds_matches)[[1]]
        
        for(cds_text in cds_texts) {
          # Extract product name
          product_match <- regexpr('/product="([^"]*)"', cds_text, perl = TRUE)
          product <- if(product_match > 0) {
            gsub('/product="([^"]*)"', "\\1", regmatches(cds_text, product_match))
          } else {
            "polyprotein"
          }
          
          # Extract protein ID
          protein_id_match <- regexpr('/protein_id="([^"]*)"', cds_text, perl = TRUE)
          protein_id <- if(protein_id_match > 0) {
            gsub('/protein_id="([^"]*)"', "\\1", regmatches(cds_text, protein_id_match))
          } else {
            NULL
          }
          
       # Only process polyproteins
          if(grepl("polyprotein", product, ignore.case = TRUE)) {
    cat("Found polyprotein:", product, "\n")
    if(is.null(protein_id)) {
      cat("WARNING: No protein_id found for this polyprotein. Skipping.\n")
    } else {
      cat("Found protein_id:", protein_id, "\n")
            # Get protein sequence using protein_id
            tryCatch({
              protein_record <- entrez_fetch(
                db = "protein", 
                id = protein_id, 
                rettype = "fasta", 
                retmode = "text"
              )
              
              # Extract sequence
              seq_lines <- strsplit(protein_record, "\n")[[1]]
              header_line <- seq_lines[1]
              sequence <- paste(seq_lines[-1], collapse = "")
              
              # Create custom header
              custom_header <- create_fasta_header(accession, taxonomy, organism, product)
              
              # Write to file
              cat(paste0(custom_header, "\n", sequence, "\n"), 
                  file = output_file, 
                  append = TRUE)
              
              polyprotein_count <- polyprotein_count + 1
              if(polyprotein_count %% 100 == 0) {
                cat("Downloaded", polyprotein_count, "polyproteins so far...\n")
              }
            }, error = function(e) {
              cat("Error fetching protein", protein_id, ":", conditionMessage(e), "\n")
            })
          }
        }
      }
    }
    # Avoid overloading the server
    Sys.sleep(1)
}
  
  cat("\nDownload complete! Added", polyprotein_count, "viral polyproteins to", output_file, "\n")
  return(output_file)
}

# If this script is run directly (not sourced), download polyproteins
if(!interactive()) {
  cat("Starting viral polyprotein download...\n")
  
  # Check for command line arguments (max_records)
  args <- commandArgs(trailingOnly = TRUE)
  max_records <- NULL
  if(length(args) > 0 && !is.na(as.numeric(args[1]))) {
    max_records <- as.numeric(args[1])
  }
  
  output_file <- download_viral_polyproteins(max_records = max_records)
  cat("Viral polyprotein download complete.\n")
  cat("Output saved to:", output_file, "\n")
}
