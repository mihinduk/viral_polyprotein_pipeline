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
  #' @return List containing output file path and stats
  #' @export
  download_viral_polyproteins <- function(credentials = NULL, max_records = NULL, output_file = NULL) {
    # Initialize result counters at the beginning
    result <- list(
      output_file = NULL,
      polyprotein_count = 0,
      skipped_count = 0
    )

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
    result$output_file <- output_file

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
      return(result)
    }

    cat("Found", total_count, "viral polyprotein records.\n")

    # Limit records if requested
    if(!is.null(max_records) && max_records < total_count) {
      cat("Limiting download to", max_records, "records.\n")
      total_count <- max_records
    }

    # Initialize counters
    polyprotein_count <- 0
    skipped_count <- 0
    batch_size <- 100

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

        # Extract definition, preserving spaces
        definition_match <- regexpr("DEFINITION\\s+([^\\n]+(?:\\n\\s+[^\\n]+)*)", record_text, perl = TRUE)
        definition <- if(definition_match > 0) {
          # Get the full multi-line definition and remove newlines and extra whitespace
          definition_text <- regmatches(record_text, definition_match)
          gsub("DEFINITION\\s+", "", definition_text)
          definition_cleaned <- gsub("\\n\\s+", " ", definition_text)
          definition_cleaned <- gsub("DEFINITION\\s+", "", definition_cleaned)
          trimws(definition_cleaned)
        } else {
          "unknown"
        }
        cat("Processing record:", accession, "-", substring(definition, 1, 50), "...\n")
        organism_match <- regexpr("ORGANISM\\s+([^\\n]+)", record_text, perl = TRUE)
        organism <- if(organism_match > 0) {
        trimws(gsub("ORGANISM\\s+", "", regmatches(record_text, organism_match)))
       
        # Extract taxonomy
        taxonomy_block <- gsub(".*ORGANISM[^\\n]*\\n\\s*([^/]*?)\\n/.*", "\\1", record_text, perl = TRUE)
        taxonomy <- gsub("\\s+", " ", taxonomy_block)

        # Try two approaches to find polyprotein CDS features
        cat("  Looking for polyprotein features in record", accession, "...\n")

        # First approach: standard regex pattern
        cds_pattern <- "/CDS\\s+([^/]+).*?/product=\"([^\"]*polyprotein[^\"]*)\".*?(/protein_id=\"([^\"]+)\")?.*?"
        cds_matches <- gregexpr(cds_pattern, record_text, perl = TRUE, ignore.case = TRUE)

        # Process CDS features from standard pattern
        polyprotein_features <- list()

        if(cds_matches[[1]][1] != -1) {
          cds_texts <- regmatches(record_text, cds_matches)[[1]]
          cat("  Found", length(cds_texts), "potential polyprotein features with standard pattern.\n")

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

            if(!is.null(protein_id)) {
              polyprotein_features[[protein_id]] <- list(product = product, protein_id = protein_id)
            }
          }
        }

        # Second approach: manually look for CDS sections and check for polyprotein product
        if(length(polyprotein_features) == 0) {
          cat("  Trying alternative approach to find polyprotein features...\n")

          # Find all CDS sections
          all_cds_pattern <- "/CDS\\s+([^/]+)"
          all_cds_matches <- gregexpr(all_cds_pattern, record_text, perl = TRUE)

          if(all_cds_matches[[1]][1] != -1) {
            # Get the positions where CDS sections start
            cds_start_positions <- attr(all_cds_matches[[1]], "capture.start")
            cds_lengths <- attr(all_cds_matches[[1]], "capture.length")

            for(i in 1:length(all_cds_matches[[1]])) {
              # Extract a chunk of text starting from the CDS
              pos <- all_cds_matches[[1]][i]
              chunk_end <- min(pos + 1000, nchar(record_text))  # Look at next 1000 chars
              chunk <- substr(record_text, pos, chunk_end)

              # Check if this chunk contains polyprotein product
              if(grepl('/product="[^"]*polyprotein[^"]*"', chunk, ignore.case = TRUE)) {
                # Extract product
                product_match <- regexpr('/product="([^"]*)"', chunk, perl = TRUE)
                product <- if(product_match > 0) {
                  gsub('/product="([^"]*)"', "\\1", regmatches(chunk, product_match))
                } else {
                  "polyprotein"
                }

                # Extract protein ID if present
                protein_id_match <- regexpr('/protein_id="([^"]*)"', chunk, perl = TRUE)
                protein_id <- if(protein_id_match > 0) {
                  gsub('/protein_id="([^"]*)"', "\\1", regmatches(chunk, protein_id_match))
                } else {
                  NULL
                }

                if(!is.null(protein_id)) {
                  cat("  Found polyprotein with alternative approach: protein_id =", protein_id, "\n")
                  polyprotein_features[[protein_id]] <- list(product = product, protein_id = protein_id)
                }
              }
            }
          }
        }

        # Process all found polyprotein features
        cat("  Total polyprotein features found:", length(polyprotein_features), "\n")

        for(feature in polyprotein_features) {
          product <- feature$product
          protein_id <- feature$protein_id

          # Only process polyproteins
          if(grepl("polyprotein", product, ignore.case = TRUE) && !is.null(protein_id)) {
                # Get protein sequence using protein_id
                tryCatch({
                  cat("  Fetching protein sequence for", protein_id, "...\n")
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
                  result$polyprotein_count <- result$polyprotein_count + 1

                  if(polyprotein_count %% 100 == 0) {
                    cat("Downloaded", polyprotein_count, "polyproteins so far...\n")
                  }
                }, error = function(e) {
                  cat("Error fetching protein", protein_id, ":", conditionMessage(e), "\n")
                  skipped_count <- skipped_count + 1
                  result$skipped_count <- result$skipped_count + 1
                })
              }
            }
          }
        }
      }
      # Avoid overloading the server
      Sys.sleep(1)
    }  # End of for loop

    cat("\nDownload complete! Added", result$polyprotein_count, "viral polyproteins to", output_file, "\n")
    cat("Skipped", result$skipped_count, "polyproteins due to missing protein IDs or fetch errors.\n")
    return(result)
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

    result <- download_viral_polyproteins(max_records = max_records)
    cat("Viral polyprotein download complete.\n")
    cat("Output saved to:", result$output_file, "\n")
    cat("Total polyproteins downloaded:", result$polyprotein_count, "\n")
  }
