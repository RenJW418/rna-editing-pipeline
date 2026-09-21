#!/usr/bin/env Rscript

usage <- function() {
  cat(paste0(
    "Read-level dinucleotide analysis from coordinate-sorted BAM files\n\n",
    "Usage:\n",
    "  Rscript dinucleotide_analysis.R --input-dir DIR --target-start INT [options]\n\n",
    "Required arguments:\n",
    "  --input-dir DIR             Directory containing BAM files\n",
    "  --target-start INT          First 1-based reference position\n\n",
    "Optional arguments:\n",
    "  --target-end INT            Last reference position (default: target-start + 1)\n",
    "  --output FILE               Output CSV (default: DIR/dinucleotide_summary.csv)\n",
    "  --bam-pattern REGEX         BAM filename pattern (default: \\.bam$)\n",
    "  --categories LIST           Comma-separated categories (default: AA,AG,GG,GA)\n",
    "  --reference NAME            Restrict analysis to one reference sequence\n",
    "  --min-mapq INT              Minimum mapping quality (default: 0)\n",
    "  --yield-size INT            Alignments read per chunk (default: 100000)\n",
    "  --include-secondary BOOL    Include secondary alignments (default: false)\n",
    "  --include-supplementary BOOL Include supplementary alignments (default: false)\n",
    "  --include-qcfail BOOL       Include vendor-QC-failed reads (default: false)\n",
    "  --exclude-duplicates BOOL   Exclude duplicate-marked reads (default: false)\n",
    "  --help                      Show this help message\n\n",
    "Boolean values: true/false, yes/no, or 1/0.\n"
  ))
}

parse_bool <- function(value, option_name) {
  normalized <- tolower(trimws(value))
  if (normalized %in% c("true", "yes", "1")) return(TRUE)
  if (normalized %in% c("false", "no", "0")) return(FALSE)
  stop(sprintf("%s must be true or false; received: %s", option_name, value), call. = FALSE)
}

parse_integer <- function(value, option_name, minimum = NULL) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || as.character(parsed) != trimws(value)) {
    stop(sprintf("%s must be an integer; received: %s", option_name, value), call. = FALSE)
  }
  if (!is.null(minimum) && parsed < minimum) {
    stop(sprintf("%s must be at least %d", option_name, minimum), call. = FALSE)
  }
  parsed
}

parse_args <- function(args) {
  config <- list(
    input_dir = NULL,
    target_start = NULL,
    target_end = NULL,
    output = NULL,
    bam_pattern = "\\.bam$",
    categories = c("AA", "AG", "GG", "GA"),
    reference = NULL,
    min_mapq = 0L,
    yield_size = 100000L,
    include_secondary = FALSE,
    include_supplementary = FALSE,
    include_qcfail = FALSE,
    exclude_duplicates = FALSE
  )

  if (length(args) == 0L || any(args %in% c("--help", "-h"))) {
    usage()
    quit(status = 0L)
  }

  known_options <- c(
    "input-dir", "target-start", "target-end", "output", "bam-pattern",
    "categories", "reference", "min-mapq", "yield-size",
    "include-secondary", "include-supplementary", "include-qcfail",
    "exclude-duplicates"
  )

  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected positional argument: %s", token), call. = FALSE)
    }

    option_text <- substring(token, 3L)
    if (grepl("=", option_text, fixed = TRUE)) {
      pieces <- strsplit(option_text, "=", fixed = TRUE)[[1L]]
      option_name <- pieces[[1L]]
      option_value <- paste(pieces[-1L], collapse = "=")
      index <- index + 1L
    } else {
      option_name <- option_text
      if (index == length(args)) {
        stop(sprintf("Missing value for --%s", option_name), call. = FALSE)
      }
      option_value <- args[[index + 1L]]
      index <- index + 2L
    }

    if (!option_name %in% known_options) {
      stop(sprintf("Unknown option: --%s", option_name), call. = FALSE)
    }

    key <- gsub("-", "_", option_name, fixed = TRUE)
    if (option_name %in% c(
      "include-secondary", "include-supplementary", "include-qcfail",
      "exclude-duplicates"
    )) {
      config[[key]] <- parse_bool(option_value, paste0("--", option_name))
    } else if (option_name %in% c("target-start", "target-end")) {
      config[[key]] <- parse_integer(option_value, paste0("--", option_name), 1L)
    } else if (option_name == "min-mapq") {
      config[[key]] <- parse_integer(option_value, "--min-mapq", 0L)
    } else if (option_name == "yield-size") {
      config[[key]] <- parse_integer(option_value, "--yield-size", 1L)
    } else if (option_name == "categories") {
      config[[key]] <- toupper(trimws(strsplit(option_value, ",", fixed = TRUE)[[1L]]))
    } else {
      config[[key]] <- option_value
    }
  }

  if (is.null(config$input_dir) || !nzchar(config$input_dir)) {
    stop("--input-dir is required", call. = FALSE)
  }
  if (is.null(config$target_start)) {
    stop("--target-start is required", call. = FALSE)
  }
  if (is.null(config$target_end)) {
    config$target_end <- config$target_start + 1L
  }
  if (config$target_end < config$target_start) {
    stop("--target-end must be greater than or equal to --target-start", call. = FALSE)
  }

  window_width <- config$target_end - config$target_start + 1L
  if (length(config$categories) == 0L || any(!nzchar(config$categories))) {
    stop("--categories must contain at least one category", call. = FALSE)
  }
  if (anyDuplicated(config$categories)) {
    stop("--categories must not contain duplicates", call. = FALSE)
  }
  if (any(nchar(config$categories) != window_width) ||
      any(!grepl("^[ACGTN]+$", config$categories))) {
    stop(sprintf(
      "Every category must contain %d A/C/G/T/N characters to match the target window",
      window_width
    ), call. = FALSE)
  }

  config$input_dir <- normalizePath(config$input_dir, mustWork = TRUE)
  if (!dir.exists(config$input_dir)) {
    stop(sprintf("Input directory does not exist: %s", config$input_dir), call. = FALSE)
  }
  if (is.null(config$output) || !nzchar(config$output)) {
    config$output <- file.path(config$input_dir, "dinucleotide_summary.csv")
  }
  config$output <- normalizePath(config$output, mustWork = FALSE)
  config
}

extract_reference_window <- function(sequence, alignment_start, cigar, target_positions) {
  if (is.na(sequence) || is.na(alignment_start) || is.na(cigar) || cigar == "*") {
    return(NA_character_)
  }

  length_matches <- regmatches(cigar, gregexpr("[0-9]+", cigar, perl = TRUE))[[1L]]
  operation_matches <- regmatches(cigar, gregexpr("[MIDNSHP=X]", cigar, perl = TRUE))[[1L]]
  if (length(length_matches) == 0L || length(length_matches) != length(operation_matches)) {
    return(NA_character_)
  }

  operation_lengths <- suppressWarnings(as.integer(length_matches))
  if (anyNA(operation_lengths)) return(NA_character_)

  query_position <- 1L
  reference_position <- as.integer(alignment_start)
  extracted <- rep(NA_character_, length(target_positions))

  for (operation_index in seq_along(operation_matches)) {
    operation <- operation_matches[[operation_index]]
    operation_length <- operation_lengths[[operation_index]]

    if (operation %in% c("M", "=", "X")) {
      covered <- which(
        target_positions >= reference_position &
          target_positions < reference_position + operation_length
      )
      if (length(covered) > 0L) {
        query_indices <- query_position + target_positions[covered] - reference_position
        extracted[covered] <- substring(sequence, query_indices, query_indices)
      }
      query_position <- query_position + operation_length
      reference_position <- reference_position + operation_length
    } else if (operation %in% c("I", "S")) {
      query_position <- query_position + operation_length
    } else if (operation %in% c("D", "N")) {
      deleted_target <- target_positions >= reference_position &
        target_positions < reference_position + operation_length
      if (any(deleted_target)) return(NA_character_)
      reference_position <- reference_position + operation_length
    } else if (!operation %in% c("H", "P")) {
      return(NA_character_)
    }

    if (all(!is.na(extracted))) break
  }

  if (anyNA(extracted) || any(nchar(extracted) != 1L)) return(NA_character_)
  paste0(extracted, collapse = "")
}

flag_is_set <- function(flags, bit) {
  bitwAnd(as.integer(flags), as.integer(bit)) != 0L
}

sample_name_from_path <- function(path) {
  name <- basename(path)
  name <- sub("\\.bam$", "", name, ignore.case = TRUE)
  sub("\\.sorted$", "", name, ignore.case = TRUE)
}

analyze_bam <- function(bam_path, config) {
  target_positions <- seq.int(config$target_start, config$target_end)
  category_counts <- setNames(rep(0, length(config$categories)), config$categories)
  total_reads <- 0
  spanning_reads <- 0

  bam <- Rsamtools::BamFile(bam_path, yieldSize = config$yield_size)
  param <- Rsamtools::ScanBamParam(
    what = c("flag", "rname", "pos", "mapq", "cigar", "seq")
  )

  opened <- FALSE
  on.exit(if (opened) close(bam), add = TRUE)
  open(bam)
  opened <- TRUE

  repeat {
    records <- Rsamtools::scanBam(bam, param = param)[[1L]]
    record_count <- length(records$pos)
    if (record_count == 0L) break

    keep <- !is.na(records$pos) & !is.na(records$mapq) & !is.na(records$flag)
    keep <- keep & !flag_is_set(records$flag, 0x4L)
    keep <- keep & records$mapq >= config$min_mapq
    if (!config$include_secondary) keep <- keep & !flag_is_set(records$flag, 0x100L)
    if (!config$include_supplementary) keep <- keep & !flag_is_set(records$flag, 0x800L)
    if (!config$include_qcfail) keep <- keep & !flag_is_set(records$flag, 0x200L)
    if (config$exclude_duplicates) keep <- keep & !flag_is_set(records$flag, 0x400L)
    if (!is.null(config$reference) && nzchar(config$reference)) {
      keep <- keep & !is.na(records$rname) & as.character(records$rname) == config$reference
    }

    selected <- which(keep)
    total_reads <- total_reads + length(selected)
    if (length(selected) == 0L) next

    for (read_index in selected) {
      dinucleotide <- extract_reference_window(
        sequence = as.character(records$seq[[read_index]]),
        alignment_start = records$pos[[read_index]],
        cigar = records$cigar[[read_index]],
        target_positions = target_positions
      )
      if (is.na(dinucleotide)) next

      dinucleotide <- toupper(dinucleotide)
      spanning_reads <- spanning_reads + 1
      if (dinucleotide %in% config$categories) {
        category_counts[[dinucleotide]] <- category_counts[[dinucleotide]] + 1
      }
    }
  }

  close(bam)
  opened <- FALSE

  filtered_reads <- sum(category_counts)
  result <- data.frame(
    sample = sample_name_from_path(bam_path),
    bam_file = normalizePath(bam_path, mustWork = TRUE),
    total_reads = total_reads,
    spanning_reads = spanning_reads,
    filtered_reads = filtered_reads,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (category in config$categories) {
    result[[paste0(category, "_count")]] <- unname(category_counts[[category]])
  }
  for (category in config$categories) {
    result[[paste0(category, "_percent")]] <- if (filtered_reads > 0) {
      100 * unname(category_counts[[category]]) / filtered_reads
    } else {
      NA_real_
    }
  }
  result
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop(
      "The Rsamtools package is required. Install it with BiocManager::install('Rsamtools').",
      call. = FALSE
    )
  }

  bam_files <- sort(list.files(
    config$input_dir,
    pattern = config$bam_pattern,
    full.names = TRUE
  ))
  if (length(bam_files) == 0L) {
    stop(sprintf(
      "No BAM files matched pattern %s in %s",
      shQuote(config$bam_pattern),
      config$input_dir
    ), call. = FALSE)
  }

  output_dir <- dirname(config$output)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop(sprintf("Could not create output directory: %s", output_dir), call. = FALSE)
  }

  message(sprintf("Found %d BAM file(s)", length(bam_files)))
  results <- lapply(seq_along(bam_files), function(index) {
    message(sprintf("[%d/%d] Processing %s", index, length(bam_files), basename(bam_files[[index]])))
    analyze_bam(bam_files[[index]], config)
  })
  summary_table <- do.call(rbind, results)
  utils::write.csv(summary_table, config$output, row.names = FALSE, na = "")
  message(sprintf("Wrote summary to %s", config$output))
}

if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(error) {
      message("ERROR: ", conditionMessage(error))
      quit(status = 1L)
    }
  )
}
