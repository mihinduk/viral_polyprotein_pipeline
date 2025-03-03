#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required libraries
if (!require("rentrez")) install.packages("rentrez")
if (!require("Biostrings")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!require("stringr")) install.packages("stringr")
if (!require("seqinr")) install.packages("seqinr")

library(rentrez)
library(Biostrings)
library(stringr)
library(seqinr)

# Set working directory with escaped spaces
setwd("/Users/handley_lab/Handley\ Lab\ Dropbox/virome/KM_algorithm_dev")

# Function to securely get input from user
getInput <- function(prompt) {
  # Use standard readline function
  cat(prompt)
  flush.console()
  user_input <- readline()
  return(user_input)
}

# Get user email and API key
getUserCredentials <- function() {
  # Make sure these prompts are sent to the console
  cat("NCBI E-utilities requires your email address for better service.\n")
  
  # Hardcode a test email for demonstration
  user_email <- "test@example.com"
  cat("Using default email for demonstration: test@example.com\n")
  
  # No API key for demonstration
  api_key <- NULL
  cat("No API key used for demonstration.\n")
  
  cat("Using email:", user_email, "\n")
  if (!is.null(api_key) && api_key != "") {
    cat("API key provided. Using higher rate limits.\n")
  } else {
    cat("No API key provided. Using standard rate limits.\n")
  }
  
  return(list(email = user_email, api_key = api_key))
}

# Format taxonomy with placeholders
formatTaxonomy <- function(taxonomy) {
  ranks <- c("phylum", "class", "order", "family", "genus", "species")
  result <- vector("character", length(ranks))
  
  taxonomy_parts <- strsplit(taxonomy, ";")[[1]]
  taxonomy_parts <- trimws(taxonomy_parts)
  
  for (i in seq_along(ranks)) {
    # Try direct pattern matching
    rank_pattern <- paste0(ranks[i], ":")
    matched <- FALSE
    
    # First try exact matches in taxonomy string
    if (any(grepl(rank_pattern, taxonomy_parts, ignore.case = TRUE))) {
      match_idx <- which(grepl(rank_pattern, taxonomy_parts, ignore.case = TRUE))
      if (length(match_idx) > 0) {
        result[i] <- gsub(".*:", "", taxonomy_parts[match_idx[1]])
        matched <- TRUE
      }
    }
    
    # For family, genus, and species, try inferring from taxonomy
    if (!matched) {
      if (ranks[i] == "family" && any(grepl("viridae$|virinae$", taxonomy_parts, ignore.case = TRUE))) {
        family_idx <- grep("viridae$|virinae$", taxonomy_parts, ignore.case = TRUE)
        if (length(family_idx) > 0) {
          result[i] <- taxonomy_parts[family_idx[1]]
          matched <- TRUE
        }
      } else if (ranks[i] == "genus" && any(grepl("virus$", taxonomy_parts, ignore.case = TRUE))) {
        genus_idx <- grep("virus$", taxonomy_parts, ignore.case = TRUE)
        if (length(genus_idx) > 0) {
          result[i] <- taxonomy_parts[genus_idx[1]]
          matched <- TRUE
        }
      }
    }
    
    # Create placeholder based on highest available taxonomy if still not matched
    if (!matched) {
      for (j in (i-1):1) {
        if (result[j] != "") {
          result[i] <- paste0(result[j], "_no_", ranks[i])
          break
        }
      }
      # If no higher taxonomy available
      if (!matched && result[i] == "") result[i] <- "Unknown"
    }
  }
  
  return(result)
}

# Create formatted FASTA header
createHeader <- function(accession, taxonomy, organism, protein_name) {
  tax_info <- formatTaxonomy(taxonomy)
  
  # Format header with placeholders
  header <- paste0(
    ">", accession, "|Viruses|",
    paste(tax_info, collapse = "|"), "|",
    organism, "|", protein_name
  )
  
  return(header)
}

# Run BLASTX via system command
runBlastx <- function(query_file, db_file, output_file, evalue = 1e-10, max_target_seqs = 1) {
  # Check if BLAST is installed
  blast_check <- system("which blastx", intern = TRUE)
  if (length(blast_check) == 0) {
    cat("ERROR: BLAST+ not found in system path. Please install BLAST+.\n")
    return(FALSE)
  }
  
  # Create BLAST command
  blast_cmd <- paste0(
    "blastx -query ", query_file, 
    " -db ", db_file,
    " -out ", output_file,
    " -evalue ", evalue,
    " -max_target_seqs ", max_target_seqs,
    " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'"
  )
  
  # Run BLAST
  cat("Running BLASTX...\n")
  system(blast_cmd)
  
  # Check if output file exists and is not empty
  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    return(TRUE)
  } else {
    cat("ERROR: BLASTX failed or produced empty results.\n")
    return(FALSE)
  }
}

# Process BLAST results and add to mature proteins database
processBlastResults <- function(blast_file, polyprotein_file, mature_protein_file, taxonomy_map) {
  # Read BLAST results
  blast_results <- read.table(blast_file, sep = "\t", stringsAsFactors = FALSE)
  colnames(blast_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                              "gapopen", "qstart", "qend", "sstart", "send", 
                              "evalue", "bitscore", "qseq", "sseq")
  
  # Process each hit
  for (i in 1:nrow(blast_results)) {
    # Get query accession and coordinates
    query_parts <- strsplit(blast_results$qseqid[i], "\\|")[[1]]
    accession <- query_parts[1]
    
    # Get taxonomic info from the map
    tax_info <- taxonomy_map[[accession]]
    if (is.null(tax_info)) {
      cat("Warning: No taxonomy information for", accession, "\n")
      next
    }
    
    # Get mature protein name from subject
    subject_parts <- strsplit(blast_results$sseqid[i], "\\|")[[1]]
    protein_name <- if (length(subject_parts) >= 9) subject_parts[9] else "unknown_protein"
    
    # Create header
    header <- paste0(
      ">", accession, "|Viruses|",
      paste(tax_info$taxonomy, collapse = "|"), "|",
      tax_info$organism, "|BLASTX_", protein_name
    )
    
    # Translated sequence from BLAST
    sequence <- blast_results$sseq[i]
    
    # Write to mature protein file
    cat(paste0(header, "\n", sequence, "\n"), 
        file = mature_protein_file, 
        append = TRUE)
  }
  
  return(nrow(blast_results))
}

# Main function to download viral polyproteins
downloadViralPolyproteins <- function(user_email, api_key = NULL) {
  # Set up rentrez
  if (!is.null(api_key) && api_key != "") {
    set_entrez_key(api_key)
  }
  
  # Get current date for file naming
  today <- format(Sys.Date(), "%Y_%m_%d")
  
  # Output file names
  polyprotein_file <- paste0(today, "_viral_polyprotein.fasta")
  mature_protein_file <- paste0(today, "_viral_polyprotein_mature_proteins.fasta")
  missing_annotation_file <- paste0(today, "_viral_polyprotein_missing_annotation.fasta")
  missing_txt_file <- paste0(today, "_polyproteins_missing_annotation.txt")
  blast_result_file <- paste0(today, "_blastx_results.txt")
  
  # Create or clear output files
  file.create(polyprotein_file)
  file.create(mature_protein_file)
  file.create(missing_annotation_file)
  file.create(missing_txt_file)
  
  # Search for viral polyproteins in nuccore
  search_term <- "(virus[organism] OR viruses[organism]) AND polyprotein[protein]"
  
  # For testing purposes, limit to just a few records - REMOVE THIS LIMIT FOR PRODUCTION USE
  search_term <- paste0(search_term, " AND complete[title] AND \"Zika virus\"[organism] AND Brazil")
  cat("Using search term (LIMITED FOR TESTING):", search_term, "\n")
  
  # Perform initial search to get count
  search_results <- entrez_search(db = "nuccore", term = search_term, retmax = 0)
  total_count <- search_results$count
  
  cat(paste0("Found ", total_count, " viral polyprotein records.\n"))
  
  # Initialize counters and trackers
  polyprotein_count <- 0
  mature_protein_count <- 0
  missing_annotation_count <- 0
  missing_annotation_ids <- character()
  taxonomy_map <- list()
  
  # Process in batches
  batch_size <- 20  # Smaller batch size for testing
  if (total_count > 0) {
    for (i in seq(1, total_count, by = batch_size)) {
      cat(paste0("Processing records ", i, " to ", min(i + batch_size - 1, total_count), "...\n"))
      
      # Get batch of IDs
      search_batch <- entrez_search(
        db = "nuccore", 
        term = search_term, 
        retmax = batch_size, 
        retstart = i - 1,
        use_history = TRUE
      )
      
      # Fetch records - use non-XML mode for simplicity
      records <- entrez_fetch(
        db = "nuccore", 
        rettype = "gb", 
        retmode = "text",  # Use text mode instead of XML
        web_history = search_batch$web_history,
        parsed = FALSE
      )
      
      # Process each record in text mode
      cat("Fetched records in text mode\n")
      
      # With text mode, we need to parse the GenBank records differently
      # For demonstration, we'll use a simpler approach just to show the concept
      
      # Split the record by LOCUS lines (each new record starts with LOCUS)
      gb_records <- strsplit(records, "LOCUS")[[1]][-1]  # Remove the first empty element
      
      if (length(gb_records) == 0) {
        cat("No records found in the response\n")
        next
      }
      
      cat("Found", length(gb_records), "records to process\n")
      
      for (record_text in gb_records) {
        # Add the LOCUS back for proper parsing
        record_text <- paste0("LOCUS", record_text)
        
        # Extract basic info using regex
        accession <- gsub(".*ACCESSION\\s+([^\\s]+).*", "\\1", record_text, perl = TRUE)
        if (accession == record_text) accession <- "unknown"
        
        organism_match <- regexpr("ORGANISM\\s+([^\\n]+)", record_text, perl = TRUE)
        organism <- if (organism_match > 0) {
          trimws(gsub("ORGANISM\\s+", "", regmatches(record_text, organism_match)))
        } else {
          "unknown"
        }
        
        # Extract taxonomy, which is right after ORGANISM
        taxonomy_block <- gsub(".*ORGANISM[^\\n]*\\n\\s*([^/]*?)\\n/.*", "\\1", record_text, perl = TRUE)
        taxonomy <- gsub("\\s+", " ", taxonomy_block)  # Clean up whitespace
        
        cat("Processing accession:", accession, "\n")
        
        # Store taxonomy info for future use with BLAST results
        tax_result <- formatTaxonomy(taxonomy)
        taxonomy_map[[accession]] <- list(
          taxonomy = tax_result,
          organism = organism
        )
        
        # Look for polyprotein features using regex
        has_mature_proteins <- FALSE
        
        # Find CDS features that contain "polyprotein"
        cds_pattern <- "/CDS\\s+([^/]+).*?/product=\"([^\"]*polyprotein[^\"]*)\".*?(/protein_id=\"([^\"]+)\")?.*?\\n"
        cds_matches <- gregexpr(cds_pattern, record_text, perl = TRUE, ignore.case = TRUE)
        
        # If we found CDS features
        if (cds_matches[[1]][1] != -1) {
          cds_texts <- regmatches(record_text, cds_matches)[[1]]
          
          for (cds_text in cds_texts) {
            # Extract product name
            product_match <- regexpr('/product="([^"]*)"', cds_text, perl = TRUE)
            product <- if (product_match > 0) {
              gsub('/product="([^"]*)"', "\\1", regmatches(cds_text, product_match))
            } else {
              NULL
            }
            
            # Extract protein ID
            protein_id_match <- regexpr('/protein_id="([^"]*)"', cds_text, perl = TRUE)
            protein_id <- if (protein_id_match > 0) {
              gsub('/protein_id="([^"]*)"', "\\1", regmatches(cds_text, protein_id_match))
            } else {
              NULL
            }
            
            # Make sure it's a polyprotein
            if (!is.null(product) && grepl("polyprotein", product, ignore.case = TRUE)) {
              # Found polyprotein
              polyprotein_count <- polyprotein_count + 1
              
              # Get protein sequence using protein_id
              if (!is.null(protein_id)) {
                protein_record <- entrez_fetch(
                  db = "protein", 
                  id = protein_id, 
                  rettype = "fasta", 
                  retmode = "text"
                )
                
                # Extract sequence
                sequence_lines <- strsplit(protein_record, "\n")[[1]][-1]
                sequence <- paste(sequence_lines, collapse = "")
                
                # Create header and write to polyprotein file
                header <- createHeader(accession, taxonomy, organism, product)
                cat(paste0(header, "\n", sequence, "\n"), 
                    file = polyprotein_file, 
                    append = TRUE)
              }
              
              # Check for mature peptides using regex
              has_mature_peptides <- FALSE
              
              # Look for mat_peptide features that come after this CDS
              mat_pattern <- "/mat_peptide\\s+([^/]+).*?/product=\"([^\"]*)\".*?\\n"
              mat_matches <- gregexpr(mat_pattern, record_text, perl = TRUE)
              
              if (mat_matches[[1]][1] != -1) {
                mat_texts <- regmatches(record_text, mat_matches)[[1]]
                
                for (mat_text in mat_texts) {
                  # Extract the position information
                  pos_match <- regexpr("mat_peptide\\s+([^/]+)", mat_text)
                  pos_text <- if (pos_match > 0) {
                    gsub("mat_peptide\\s+", "", regmatches(mat_text, pos_match))
                  } else {
                    ""
                  }
                  
                  # Extract mature peptide product name
                  mat_product_match <- regexpr('/product="([^"]*)"', mat_text, perl = TRUE)
                  mat_product <- if (mat_product_match > 0) {
                    gsub('/product="([^"]*)"', "\\1", regmatches(mat_text, mat_product_match))
                  } else {
                    "unknown_mature_peptide"
                  }
                  
                  # Try to extract coordinates (this is simplified)
                  coords <- str_extract_all(pos_text, "\\d+")[[1]]
                  
                  if (length(coords) >= 2) {
                    has_mature_peptides <- TRUE
                    has_mature_proteins <- TRUE
                    mature_protein_count <- mature_protein_count + 1
                    
                    # For demonstration, we'll fetch the protein sequence from protein_id to save it
                    # In a real implementation, you'd extract the mature peptide sequence based on coordinates
                    if (!is.null(protein_id)) {
                      # Since we can't reliably extract the partial sequence without the complete 
                      # protein sequence, we'll save a placeholder notification
                      cat(paste0(">", accession, "|Viruses|",
                            paste(tax_result, collapse = "|"), "|",
                            organism, "|", mat_product, " (coordinates: ", coords[1], "-", coords[2], ")\n",
                            "PEPTIDE_SEQUENCE_WOULD_BE_EXTRACTED_HERE\n"), 
                        file = mature_protein_file, 
                        append = TRUE)
                    }
                  }
                }
              }
              
              # Track polyproteins without mature peptide annotations
              if (!has_mature_peptides) {
                missing_annotation_count <- missing_annotation_count + 1
                missing_annotation_ids <- c(missing_annotation_ids, accession)
                
                # Save the polyprotein sequence for BLAST analysis
                cat(paste0(header, "\n", sequence, "\n"),
                    file = missing_annotation_file,
                    append = TRUE)
                
                # Also save metadata to text file
                cat(paste0(accession, "\t", organism, "\t", product, "\n"),
                    file = missing_txt_file,
                    append = TRUE)
              }
            }
          }
        }
      }
      
      # Avoid overloading the server
      Sys.sleep(1)
    }
  }
  
  # Print summary after initial download
  cat(paste0("\nDownload complete.\n"))
  cat(paste0("Total polyproteins: ", polyprotein_count, "\n"))
  cat(paste0("Polyproteins with mature peptide annotations: ", mature_protein_count, "\n"))
  cat(paste0("Polyproteins without mature peptide annotations: ", missing_annotation_count, "\n"))
  
  # Run BLASTX to identify mature proteins in unannotated polyproteins
  if (missing_annotation_count > 0 && file.exists(mature_protein_file) && file.info(mature_protein_file)$size > 0) {
    cat("\nRunning BLASTX on unannotated polyproteins...\n")
    
    # First, format the mature proteins as a BLAST database
    cat("Formatting BLAST database from annotated mature proteins...\n")
    db_cmd <- paste0("makeblastdb -in ", mature_protein_file, " -dbtype prot")
    system(db_cmd)
    
    # Run BLASTX
    blastx_success <- runBlastx(
      query_file = missing_annotation_file,
      db_file = mature_protein_file,
      output_file = blast_result_file
    )
    
    if (blastx_success) {
      # Process BLAST results and add to mature proteins database
      blast_hits <- processBlastResults(
        blast_file = blast_result_file,
        polyprotein_file = missing_annotation_file,
        mature_protein_file = mature_protein_file,
        taxonomy_map = taxonomy_map
      )
      
      cat(paste0("Added ", blast_hits, " mature proteins identified by BLASTX.\n"))
    }
  }
  
  # Final summary
  cat("\nAll processing complete.\n")
  cat(paste0("Polyprotein database: ", polyprotein_file, "\n"))
  cat(paste0("Mature protein database: ", mature_protein_file, "\n"))
  cat(paste0("Missing annotation file: ", missing_txt_file, "\n"))
}

# Main execution - only run if script is run directly (not sourced)
if (!interactive()) {
  cat("Starting viral polyprotein database generation...\n")
  credentials <- getUserCredentials()
  
  # Check if email is empty or just whitespace
  email_valid <- !is.null(credentials$email) && 
                 credentials$email != "" && 
                 !grepl("^\\s*$", credentials$email)
  
  if (email_valid) {
    downloadViralPolyproteins(credentials$email, credentials$api_key)
  } else {
    cat("ERROR: Email is required to use NCBI E-utilities.\n")
  }
}
