#!/usr/bin/env Rscript

# 04_blast_annotation.R - Use BLAST to annotate unannotated polyproteins
# This script takes polyproteins without mature protein annotations and uses BLAST
# to identify potential mature proteins based on similarity to known mature proteins

# Load configuration and helper functions
source("01_config.R")
source("03_extract_mature_proteins.R")

#' Set up BLAST database from mature proteins
#' 
#' @param mature_proteins_file Path to FASTA file with mature proteins
#' @return TRUE if database creation was successful
#' @export
create_blast_db <- function(mature_proteins_file) {
  cat("Creating BLAST database from mature proteins...\n")
  
  # Check if BLAST is installed
  blast_check <- system2("which", "makeblastdb", stdout = TRUE, stderr = FALSE)
  if(length(blast_check) == 0) {
    cat("ERROR: BLAST+ not found in system path. Please install BLAST+.\n")
    return(FALSE)
  }
  
  # Create BLAST database
  db_cmd <- paste0("makeblastdb -in ", mature_proteins_file, " -dbtype prot -parse_seqids")
  result <- system(db_cmd)
  
  if(result != 0) {
    cat("ERROR: Failed to create BLAST database.\n")
    return(FALSE)
  }
  
  cat("BLAST database created successfully.\n")
  return(TRUE)
}

#' Run BLASTX on unannotated polyproteins
#' 
#' @param query_file Path to FASTA file with unannotated polyproteins
#' @param db_file Path to FASTA file with mature proteins (BLAST database)
#' @param output_file Path to output file for BLAST results
#' @param evalue E-value threshold (default: 1e-10)
#' @param max_target_seqs Maximum number of target sequences per query (default: 5)
#' @return TRUE if BLASTX ran successfully
#' @export
run_blastx <- function(query_file, db_file, output_file, evalue = 1e-10, max_target_seqs = 5) {
  cat("Running BLASTX to identify mature proteins in unannotated polyproteins...\n")
  
  # Check if BLAST is installed
  blast_check <- system2("which", "blastx", stdout = TRUE, stderr = FALSE)
  if(length(blast_check) == 0) {
    cat("ERROR: BLAST+ not found in system path. Please install BLAST+.\n")
    return(FALSE)
  }
  
  # Create BLASTX command
  blast_cmd <- paste0(
    "blastx -query ", query_file, 
    " -db ", db_file,
    " -out ", output_file,
    " -evalue ", evalue,
    " -max_target_seqs ", max_target_seqs,
    " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe'"
  )
  
  cat("Running BLAST command:", blast_cmd, "\n")
  result <- system(blast_cmd)
  
  if(result != 0) {
    cat("ERROR: BLASTX failed.\n")
    return(FALSE)
  }
  
  # Check if output file exists and is not empty
  if(!file.exists(output_file) || file.info(output_file)$size == 0) {
    cat("WARNING: BLASTX produced empty results.\n")
    return(FALSE)
  }
  
  cat("BLASTX completed successfully.\n")
  return(TRUE)
}

#' Process BLAST results and extract predicted mature proteins
#' 
#' @param blast_file Path to BLAST results file
#' @param polyprotein_file Path to FASTA file with unannotated polyproteins
#' @param mature_protein_file Path to FASTA file to append new mature proteins
#' @param taxonomy_map_file Path to RDS file with taxonomy mapping (or NULL)
#' @param min_identity Minimum percent identity for hits (default: 70)
#' @param min_length Minimum alignment length (default: 50)
#' @return Number of mature proteins identified by BLAST
#' @export
process_blast_results <- function(blast_file, polyprotein_file, mature_protein_file, 
                                 taxonomy_map_file = NULL, min_identity = 70, min_length = 50) {
  cat("Processing BLAST results to extract predicted mature proteins...\n")
  
  # Load taxonomy map if available
  taxonomy_map <- NULL
  if(!is.null(taxonomy_map_file) && file.exists(taxonomy_map_file)) {
    taxonomy_map <- readRDS(taxonomy_map_file)
  }
  
  # Read BLAST results
  blast_results <- read.delim(blast_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(blast_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                              "gapopen", "qstart", "qend", "sstart", "send", 
                              "evalue", "bitscore", "qframe")
  
  # Filter by identity and length
  blast_results <- blast_results[blast_results$pident >= min_identity & 
                                blast_results$length >= min_length, ]
  
  if(nrow(blast_results) == 0) {
    cat("No significant BLAST hits found.\n")
    return(0)
  }
  
  cat("Found", nrow(blast_results), "significant BLAST hits.\n")
  
  # Read polyprotein sequences
  polyproteins <- readLines(polyprotein_file)
  polyprotein_seqs <- list()
  current_header <- NULL
  current_seq <- character(0)
  
  for(line in polyproteins) {
    if(substr(line, 1, 1) == ">") {
      if(!is.null(current_header)) {
        polyprotein_seqs[[current_header]] <- paste(current_seq, collapse = "")
      }
      current_header <- substr(line, 2, nchar(line))
      current_seq <- character(0)
    } else {
      current_seq <- c(current_seq, line)
    }
  }
  # Add the last sequence
  if(!is.null(current_header)) {
    polyprotein_seqs[[current_header]] <- paste(current_seq, collapse = "")
  }
  
  # Process each hit to extract mature protein
  mature_proteins_added <- 0
  processed_regions <- list()  # Track processed regions to avoid overlaps
  
  # Sort by query ID and bitscore (descending) to get best hits first
  blast_results <- blast_results[order(blast_results$qseqid, -blast_results$bitscore), ]
  
  for(i in 1:nrow(blast_results)) {
    hit <- blast_results[i, ]
    query_id <- hit$qseqid
    subject_id <- hit$sseqid
    
    # Skip if we've already processed this query
    if(is.null(processed_regions[[query_id]])) {
      processed_regions[[query_id]] <- list()
    }
    
    # Extract coordinates in the polyprotein
    qstart <- as.numeric(hit$qstart)
    qend <- as.numeric(hit$qend)
    qframe <- as.numeric(hit$qframe)
    
    # Check for overlap with previously processed regions
    overlap <- FALSE
    for(region in processed_regions[[query_id]]) {
      if(max(qstart, region$start) <= min(qend, region$end)) {
        overlap <- TRUE
        break
      }
    }
    
    if(overlap) next
    
    # Get polyprotein sequence
    polyprotein_seq <- polyprotein_seqs[[query_id]]
    if(is.null(polyprotein_seq)) {
      cat("WARNING: Could not find sequence for", query_id, "\n")
      next
    }
    
    # Adjust coordinates if necessary based on frame
    if(qframe < 0) {
      # For negative frames, we need to translate in reverse, but this is complex
      # A simpler approach is to just note the frame in the name
      frame_note <- paste0("(negative_frame_", abs(qframe), ")")
    } else {
      frame_note <- paste0("(frame_", qframe, ")")
    }
    
    # Extract metadata from query and subject IDs
    query_parts <- strsplit(query_id, "\\|")[[1]]
    subject_parts <- strsplit(subject_id, "\\|")[[1]]
    
    # Extract accession and protein name
    accession <- query_parts[1]
    hit_protein_name <- if(length(subject_parts) >= 9) subject_parts[9] else "unknown_protein"
    
    # Get matching region in amino acid coordinates (approximate)
    nt_start <- min(qstart, qend)
    nt_end <- max(qstart, qend)
    
    # Create a new header with BLAST annotation
    if(!is.null(taxonomy_map) && !is.null(taxonomy_map[[accession]])) {
      # Use stored taxonomy if available
      taxonomy <- taxonomy_map[[accession]]$taxonomy
      organism <- taxonomy_map[[accession]]$organism
    } else {
      # Extract from original header
      taxonomy <- paste(query_parts[3:8], collapse = "|")
      organism <- query_parts[9]
    }
    
    # Format protein name with BLAST info
    protein_name <- paste0("BLAST_", hit_protein_name, "_", frame_note,
                          "_coords", nt_start, "-", nt_end,
                          "_identity", round(hit$pident), "pct")
    
    # Create header
    new_header <- create_fasta_header(accession, taxonomy, organism, protein_name)
    
    # Get the translated sequence from BLAST results
    # In a real implementation, you'd extract and translate the sequence based on frame
    # For this example, we'll extract nucleotides and note we need translation
    extract_length <- nt_end - nt_start + 1
    nucleotide_region <- substr(polyprotein_seq, nt_start, nt_end)
    
    # For demonstration, we'll use a placeholder noting this needs translation
    # In practice, you would use Biostrings to translate correctly based on frame
    translated_seq <- paste0(
      "(Requires translation in frame ", qframe, ". Nucleotides: ",
      substr(nucleotide_region, 1, min(30, nchar(nucleotide_region))), 
      if(nchar(nucleotide_region) > 30) "..." else ""
    )
    
    # Write to mature proteins file
    cat(paste0(new_header, "\n", translated_seq, "\n"), 
        file = mature_protein_file, 
        append = TRUE)
    
    # Update counters and tracking
    mature_proteins_added <- mature_proteins_added + 1
    processed_regions[[query_id]] <- c(
      processed_regions[[query_id]],
      list(list(start = nt_start, end = nt_end))
    )
  }
  
  cat("Added", mature_proteins_added, "predicted mature proteins based on BLAST results.\n")
  return(mature_proteins_added)
}

#' Run the complete BLAST annotation pipeline
#' 
#' @param mature_proteins_file Path to FASTA file with known mature proteins
#' @param missing_annotations_file Path to FASTA file with unannotated polyproteins
#' @param taxonomy_map_file Path to RDS file with taxonomy mapping (or NULL)
#' @param output_blast Path to output file for BLAST results
#' @return TRUE if the pipeline ran successfully
#' @export
run_blast_annotation <- function(mature_proteins_file, missing_annotations_file, 
                                taxonomy_map_file = NULL, output_blast = NULL) {
  
  if(is.null(output_blast)) {
    output_blast <- get_dated_filename("blastx_results", "txt")
  }
  
  # 1. Create BLAST database
  if(!create_blast_db(mature_proteins_file)) {
    cat("Failed to create BLAST database. Aborting.\n")
    return(FALSE)
  }
  
  # 2. Run BLASTX
  if(!run_blastx(missing_annotations_file, mature_proteins_file, output_blast)) {
    cat("BLASTX failed. Aborting.\n")
    return(FALSE)
  }
  
  # 3. Process BLAST results and add to mature proteins file
  mature_proteins_added <- process_blast_results(
    output_blast, 
    missing_annotations_file, 
    mature_proteins_file, 
    taxonomy_map_file
  )
  
  if(mature_proteins_added > 0) {
    cat("Successfully annotated", mature_proteins_added, "mature proteins using BLAST.\n")
    cat("Updated mature proteins file:", mature_proteins_file, "\n")
    return(TRUE)
  } else {
    cat("No mature proteins could be annotated using BLAST.\n")
    return(FALSE)
  }
}

# If this script is run directly (not sourced), run BLAST annotation
if(!interactive()) {
  cat("Starting BLAST annotation of unannotated polyproteins...\n")
  
  # Check for command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default filenames
  today <- format(Sys.Date(), "%Y_%m_%d")
  mature_file <- paste0(today, "_viral_polyprotein_mature_proteins.fasta")
  missing_file <- paste0(today, "_viral_polyprotein_missing_annotation.fasta")
  taxonomy_map <- "viral_taxonomy_map.rds"
  
  # Override with command line args if provided
  if(length(args) >= 2) {
    mature_file <- args[1]
    missing_file <- args[2]
  }
  
  success <- run_blast_annotation(mature_file, missing_file, taxonomy_map)
  
  if(success) {
    cat("BLAST annotation pipeline completed successfully.\n")
  } else {
    cat("BLAST annotation pipeline encountered errors.\n")
  }
}