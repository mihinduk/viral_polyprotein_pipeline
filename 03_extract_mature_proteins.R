#!/usr/bin/env Rscript

# 03_extract_mature_proteins.R - Extract mature proteins from viral polyproteins
# This script extracts mature proteins from the GenBank records and saves them to a FASTA file

# Load configuration and helper functions
source("01_config.R")
source("02_download_polyproteins.R")

#' Extract mature proteins from viral polyprotein records
#' 
#' @param credentials NCBI credentials (email and API key)
#' @param max_records Maximum number of records to process (NULL for all)
#' @param output_file Output FASTA filename for mature proteins (NULL for auto-generated)
#' @param missing_file Output FASTA filename for polyproteins without mature protein annotations (NULL for auto-generated)
#' @return List with paths to output files
#' @export
extract_mature_proteins <- function(credentials = NULL, max_records = NULL, 
                                   output_file = NULL, missing_file = NULL) {
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
  
  # Set output filenames
  if(is.null(output_file)) {
    output_file <- get_dated_filename("viral_polyprotein_mature_proteins")
  }
  
  if(is.null(missing_file)) {
    missing_file <- get_dated_filename("viral_polyprotein_missing_annotation")
  }
  
  # Create a file to track missing annotation metadata
  missing_meta_file <- gsub("\\.fasta$", ".txt", missing_file)
  
  # Create or clear output files
  file.create(output_file)
  file.create(missing_file)
  file.create(missing_meta_file)
  
  cat("Mature proteins will be saved to:", output_file, "\n")
  cat("Polyproteins without mature annotations will be saved to:", missing_file, "\n")
  cat("Missing annotation metadata will be saved to:", missing_meta_file, "\n")
  
  # Write header for missing metadata file
  cat(paste("accession", "organism", "product", "has_protein_id", "protein_id", sep="\t"), 
      file=missing_meta_file, append=TRUE)
  cat("\n", file=missing_meta_file, append=TRUE)
  
  # Search term for viral polyproteins
  search_term <- "(virus[organism] OR viruses[organism]) AND polyprotein[protein] AND complete[title]"
  
  # Perform initial search to get count
  cat("Searching NCBI for viral polyproteins...\n")
  search_results <- entrez_search(db = "nuccore", term = search_term, retmax = 0)
  total_count <- search_results$count
  
  if(total_count == 0) {
    cat("No records found matching the search criteria.\n")
    return(list(mature_proteins = output_file, missing_annotations = missing_file))
  }
  
  cat("Found", total_count, "viral polyprotein records.\n")
  
  # Limit records if requested
  if(!is.null(max_records) && max_records < total_count) {
    cat("Limiting processing to", max_records, "records.\n")
    total_count <- max_records
  }
  
  # Initialize counters
  polyprotein_count <- 0
  mature_protein_count <- 0
  missing_annotation_count <- 0
  batch_size <- 50  # Smaller batch size for more detailed processing
  
  # Store taxonomy info for future use
  taxonomy_map <- list()
  
  # Process in batches
  for(i in seq(1, total_count, by = batch_size)) {
    batch_end <- min(i + batch_size - 1, total_count)
    cat(paste0("Processing records ", i, " to ", batch_end, " of ", total_count, "...\n"))
    
    # Get batch of IDs
    search_batch <- entrez_search(
      db = "nuccore", 
      term = search_term, 
      retmax = batch_size, 
      retstart = i - 1,
      use_history = TRUE
    )
    
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
      
      # Store for future use
      taxonomy_map[[accession]] <- list(
        taxonomy = taxonomy,
        organism = organism
      )
      
      # Find CDS features with polyproteins
      cds_pattern <- "/CDS\\s+([^/]+).*?/product=\"([^\"]*polyprotein[^\"]*)\".*?(/protein_id=\"([^\"]+)\")?.*?\\n"
      cds_matches <- gregexpr(cds_pattern, record_text, perl = TRUE, ignore.case = TRUE)
      
      # Process polyprotein features
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
            polyprotein_count <- polyprotein_count + 1
            polyprotein_sequence <- NULL
            
            # Get protein sequence if protein_id is available
            if(!is.null(protein_id)) {
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
                polyprotein_sequence <- paste(seq_lines[-1], collapse = "")
                
              }, error = function(e) {
                cat("Error fetching protein", protein_id, ":", conditionMessage(e), "\n")
              })
            }
            
            # Look for mature peptide annotations
            has_mature_peptides <- FALSE
            
            # Extract mature peptide features (mat_peptide)
            mat_pattern <- "/mat_peptide\\s+([^/]+).*?/product=\"([^\"]*)\".*?\\n"
            mat_matches <- gregexpr(mat_pattern, record_text, perl = TRUE)
            
            if(mat_matches[[1]][1] != -1 && !is.null(polyprotein_sequence)) {
              mat_texts <- regmatches(record_text, mat_matches)[[1]]
              
              for(mat_text in mat_texts) {
                # Extract position
                pos_match <- regexpr("mat_peptide\\s+([^/]+)", mat_text)
                pos_text <- if(pos_match > 0) {
                  gsub("mat_peptide\\s+", "", regmatches(mat_text, pos_match))
                } else {
                  NULL
                }
                
                # Extract product name
                mat_product_match <- regexpr('/product="([^"]*)"', mat_text, perl = TRUE)
                mat_product <- if(mat_product_match > 0) {
                  gsub('/product="([^"]*)"', "\\1", regmatches(mat_text, mat_product_match))
                } else {
                  "unknown_mature_peptide"
                }
                
                # Calculate coordinates
                if(!is.null(pos_text)) {
                  # Extract coordinates - this is a simplified approach
                  # In practice, you'd need more robust parsing for join() features
                  coords <- as.numeric(str_extract_all(pos_text, "\\d+")[[1]])
                  
                  if(length(coords) >= 2) {
                    has_mature_peptides <- TRUE
                    
                    # Calculate protein coordinates (approximate conversion from nucleotide)
                    start_pos <- ceiling(coords[1] / 3)
                    end_pos <- floor(coords[2] / 3)
                    
                    # Make sure coordinates are within range
                    if(start_pos > 0 && end_pos <= nchar(polyprotein_sequence) && start_pos <= end_pos) {
                      # Extract mature peptide sequence from polyprotein
                      mat_sequence <- substr(polyprotein_sequence, start_pos, end_pos)
                      
                      # Create custom header
                      mat_header <- create_fasta_header(accession, taxonomy, organism, mat_product)
                      
                      # Write to mature proteins file
                      cat(paste0(mat_header, "\n", mat_sequence, "\n"), 
                          file = output_file, 
                          append = TRUE)
                      
                      mature_protein_count <- mature_protein_count + 1
                    }
                  }
                }
              }
            }
            
            # Track polyproteins without mature annotations
            if(!has_mature_peptides && !is.null(polyprotein_sequence)) {
              missing_annotation_count <- missing_annotation_count + 1
              
              # Create header for missing annotation file
              missing_header <- create_fasta_header(accession, taxonomy, organism, product)
              
              # Save polyprotein to missing annotations file
              cat(paste0(missing_header, "\n", polyprotein_sequence, "\n"),
                  file = missing_file,
                  append = TRUE)
              
              # Save metadata
              cat(paste(accession, organism, product, !is.null(protein_id), 
                        ifelse(is.null(protein_id), "NA", protein_id), sep="\t"),
                  file = missing_meta_file, append = TRUE)
              cat("\n", file = missing_meta_file, append = TRUE)
            }
          }
        }
      }
    }
    
    # Avoid overloading the server
    Sys.sleep(1)
    
    # Status update
    cat("Progress:", round(batch_end/total_count*100), "% complete\n")
    cat("Polyproteins found:", polyprotein_count, "\n")
    cat("Mature proteins extracted:", mature_protein_count, "\n")
    cat("Polyproteins without mature annotations:", missing_annotation_count, "\n")
  }
  
  # Save taxonomy map for use in BLAST processing
  saveRDS(taxonomy_map, "viral_taxonomy_map.rds")
  
  cat("\nExtraction complete!\n")
  cat("Total polyproteins processed:", polyprotein_count, "\n")
  cat("Total mature proteins extracted:", mature_protein_count, "\n")
  cat("Total polyproteins without mature annotations:", missing_annotation_count, "\n")
  
  return(list(
    mature_proteins = output_file, 
    missing_annotations = missing_file,
    missing_metadata = missing_meta_file,
    taxonomy_map = "viral_taxonomy_map.rds"
  ))
}

# If this script is run directly (not sourced), extract mature proteins
if(!interactive()) {
  cat("Starting extraction of mature proteins from viral polyproteins...\n")
  
  # Check for command line arguments (max_records)
  args <- commandArgs(trailingOnly = TRUE)
  max_records <- NULL
  if(length(args) > 0 && !is.na(as.numeric(args[1]))) {
    max_records <- as.numeric(args[1])
  }
  
  output_files <- extract_mature_proteins(max_records = max_records)
  cat("Mature protein extraction complete.\n")
  cat("Mature proteins saved to:", output_files$mature_proteins, "\n")
  cat("Polyproteins without mature annotations saved to:", output_files$missing_annotations, "\n")
}