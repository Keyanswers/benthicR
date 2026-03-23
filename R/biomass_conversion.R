#' Convert Biomass Units Using Taxon-Specific Conversion Factors
#'
#' Converts biomass measurements between different units (e.g., wet mass to dry
#' mass, AFDM to energy) using taxon-specific conversion factors. The function
#' includes a comprehensive internal database of conversion factors for marine
#' benthic taxa compiled from literature.
#'
#' @param DF Data frame with biomass data. Must contain:
#'   \itemize{
#'     \item A taxonomic column (matching taxon_col)
#'     \item One or more biomass columns to convert
#'   }
#' @param taxon_col Character or numeric. Column name or number containing
#'   taxonomic group (e.g., "Phylum", "Class", "Order"). Must match column
#'   names in the conversion factors table (default = "Group").
#' @param biomass_cols Numeric vector or character vector. Columns containing
#'   biomass values to convert. Can be column numbers (e.g., 3:5) or names
#'   (e.g., c("St1", "St2")).
#' @param conversion Character. Type of conversion to apply. Options:
#'   \itemize{
#'     \item \code{"WM_to_DM"} - Wet mass to dry mass
#'     \item \code{"DM_to_AFDM"} - Dry mass to ash-free dry mass
#'     \item \code{"WM_to_AFDM"} - Wet mass to ash-free dry mass (direct)
#'     \item \code{"DM_to_Energy"} - Dry mass to energy (J/mg)
#'     \item \code{"AFDM_to_Energy"} - AFDM to energy (J/mg)
#'   }
#' @param factors Data frame with conversion factors. If NULL, uses internal
#'   database (default = NULL). If provided, must contain columns matching
#'   taxon_col and the conversion type.
#' @param na_action Character. How to handle taxa without conversion factors:
#'   \itemize{
#'     \item \code{"warning"} - Issue warning and keep original values (default)
#'     \item \code{"remove"} - Remove rows with missing factors
#'     \item \code{"mean"} - Use mean factor across all taxa
#'   }
#'
#' @return Data frame with converted biomass values. Original data frame structure
#'   is maintained, with biomass columns containing converted values.
#'
#' @details
#' \strong{Conversion Types:}
#'
#' \itemize{
#'   \item \strong{WM_to_DM:} Wet mass → Dry mass (after drying at 60°C)
#'   \item \strong{DM_to_AFDM:} Dry mass → Ash-free dry mass (after combustion)
#'   \item \strong{WM_to_AFDM:} Wet mass → AFDM (combined conversion)
#'   \item \strong{DM_to_Energy:} Dry mass → Energy content (Joules per mg DM)
#'   \item \strong{AFDM_to_Energy:} AFDM → Energy content (Joules per mg AFDM)
#' }
#'
#' \strong{Internal Database:}
#'
#' The function includes conversion factors for major benthic taxa compiled from
#' literature sources including:
#' \itemize{
#'   \item Echinodermata (Asteroidea, Echinoidea, Ophiuroidea, Holothuroidea, Crinoidea)
#'   \item Crustacea (Amphipoda, Isopoda, Decapoda, Mysidacea, Cumacea, etc.)
#'   \item Mollusca (Bivalvia, Gastropoda, Cephalopoda)
#'   \item Annelida (Polychaeta Errantia, Polychaeta Sedentaria, Oligochaeta)
#'   \item Other groups (Porifera, Cnidaria, Bryozoa, Nemertea, etc.)
#'   \item Fish (Demersal, Pelagic, Mesopelagic)
#'   \item Insect larvae (Chironomidae, Ephemeroptera, Odonata, Trichoptera)
#' }
#'
#' For taxa not in the database, the function can use mean values or keep
#' original measurements.
#'
#' \strong{Shell-Free Dry Mass (SFDM):}
#'
#' For molluscs, factors are based on shell-free dry mass. If your data includes
#' shells, pre-process to remove shell weight or adjust factors accordingly.
#'
#' @note
#' Conversion factors are approximations and can vary with season, location,
#' and individual size. For critical applications, consider using taxon-specific
#' or site-specific factors.
#'
#' @examples
#' \dontrun{
#' # Example 1: Convert wet mass to dry mass
#' biomass_data <- data.frame(
#'   Group = c("Amphipoda", "Bivalvia", "Polych. Errantia", "Asteroidea"),
#'   Species = c("Ampelisca", "Nucula", "Glycera", "Astropecten"),
#'   St1_WM = c(10.5, 25.3, 15.8, 100.2),
#'   St2_WM = c(12.3, 28.1, 14.2, 95.7),
#'   St3_WM = c(11.7, 26.8, 16.5, 102.3)
#' )
#'
#' # Check available factors first
#' view_conversion_factors(conversion = "WM_to_DM")
#'
#' biomass_DM <- ConvertBiomass(
#'   DF = biomass_data,
#'   taxon_col = "Group",
#'   biomass_cols = 3:5,  # Columns St1_WM, St2_WM, St3_WM
#'   conversion = "WM_to_DM"
#' )
#'
#' print(biomass_DM)
#'
#' # Example 2: Convert dry mass to energy content
#' # Note: Column names changed from St1_WM to St1_DM after conversion
#' colnames(biomass_DM)[3:5] <- c("St1_DM", "St2_DM", "St3_DM")
#'
#' biomass_energy <- ConvertBiomass(
#'   DF = biomass_DM,
#'   taxon_col = "Group",
#'   biomass_cols = 3:5,  # Now contains DM values
#'   conversion = "DM_to_Energy"
#' )
#'
#' # Example 3: Use custom conversion factors
#' my_factors <- data.frame(
#'   Group = c("Amphipoda", "Bivalvia"),
#'   WM_to_DM = c(0.22, 0.09),
#'   DM_to_AFDM = c(0.75, 0.85)
#' )
#'
#' biomass_custom <- ConvertBiomass(
#'   DF = biomass_data,
#'   taxon_col = "Group",
#'   biomass_cols = 3:5,
#'   conversion = "WM_to_DM",
#'   factors = my_factors
#' )
#'
#' # Example 4: Search for available groups
#' # Find all Crustacea groups
#' view_conversion_factors(group_pattern = "acea")
#'
#' # Find all Polychaeta groups
#' view_conversion_factors(group_pattern = "Polych")
#'
#' # Example 5: Handle missing taxa with mean factor
#' biomass_with_na <- data.frame(
#'   Group = c("Amphipoda", "UnknownGroup", "Bivalvia"),
#'   Species = c("Ampelisca", "Mystery", "Nucula"),
#'   St1 = c(10, 15, 20)
#' )
#'
#' # Use mean factor for unknown groups
#' result <- ConvertBiomass(
#'   DF = biomass_with_na,
#'   taxon_col = "Group",
#'   biomass_cols = 3,
#'   conversion = "WM_to_DM",
#'   na_action = "mean"
#' )
#' }
#'
#' @references
#' Brey, T., et al. (1988). Growth and production of Macoma balthica in the
#' Wadden Sea. Marine Ecology Progress Series, 45, 119-130.
#'
#' Ricciardi, A., & Bourget, E. (1998). Weight-to-weight conversion factors
#' for marine benthic macroinvertebrates. Marine Ecology Progress Series, 163, 245-251.
#'
#' @seealso \code{\link{Functional}} for functional trait analysis,
#'   \code{\link{ExploreComm}} for community composition
#'
#' @export
ConvertBiomass <- function(DF,
                           taxon_col = "Group",
                           biomass_cols,
                           conversion,
                           factors = NULL,
                           na_action = "warning") {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  if (!is.data.frame(DF)) {
    stop("DF must be a data frame")
  }

  # Get taxon column name
  if (is.numeric(taxon_col)) {
    taxon_name <- colnames(DF)[taxon_col]
  } else {
    taxon_name <- taxon_col
  }

  if (!taxon_name %in% colnames(DF)) {
    stop("Taxon column '", taxon_name, "' not found in DF")
  }

  if (missing(biomass_cols)) {
    stop("biomass_cols must be specified")
  }

  if (missing(conversion)) {
    stop("conversion type must be specified")
  }

  # Convert biomass_cols to numbers if names
  if (is.character(biomass_cols)) {
    biomass_cols_nums <- which(colnames(DF) %in% biomass_cols)
    if (length(biomass_cols_nums) == 0) {
      stop("No matching biomass columns found in DF")
    }
    biomass_cols <- biomass_cols_nums
  }

  # Validate conversion type
  valid_conversions <- c("WM_to_DM", "DM_to_AFDM", "WM_to_AFDM",
                         "DM_to_Energy", "AFDM_to_Energy")
  if (!conversion %in% valid_conversions) {
    stop("Invalid conversion type. Must be one of: ",
         paste(valid_conversions, collapse = ", "))
  }

  # Validate na_action
  if (!na_action %in% c("warning", "remove", "mean")) {
    stop("na_action must be 'warning', 'remove', or 'mean'")
  }

  # ============================================================================
  # LOAD OR USE PROVIDED FACTORS
  # ============================================================================

  if (is.null(factors)) {
    # Use internal database
    factors <- get_conversion_factors()
  } else {
    if (!is.data.frame(factors)) {
      stop("factors must be a data frame")
    }
  }

  # Map conversion type to column name
  factor_col <- switch(conversion,
                       "WM_to_DM" = "WM_to_DM",
                       "DM_to_AFDM" = "DM_to_AFDM",
                       "WM_to_AFDM" = "WM_to_AFDM",
                       "DM_to_Energy" = "J_per_mgDM",
                       "AFDM_to_Energy" = "J_per_mgAFDM"
  )

  # Validate factor column exists
  if (!factor_col %in% colnames(factors)) {
    stop("Conversion factor '", factor_col, "' not found in factors table")
  }

  # ============================================================================
  # MERGE AND CONVERT
  # ============================================================================

  # Create copy of DF to preserve original order
  DF_converted <- DF

  # Add row index to preserve order
  DF$row_index <- seq_len(nrow(DF))

  # Merge with factors
  DF_merged <- merge(DF, factors[, c("Group", factor_col)],
                     by.x = taxon_name,
                     by.y = "Group",
                     all.x = TRUE)

  # Restore original order using row_index
  DF_merged <- DF_merged[order(DF_merged$row_index), ]

  # Remove row_index column
  DF_merged$row_index <- NULL
  DF_converted$row_index <- NULL

  # Check for missing factors
  missing_taxa <- unique(DF[[taxon_name]][is.na(DF_merged[[factor_col]])])

  if (length(missing_taxa) > 0) {
    if (na_action == "warning") {
      warning("No conversion factors found for: ",
              paste(missing_taxa, collapse = ", "),
              ". Original values retained.")

      # Keep original values for taxa without factors
      # (do nothing, conversion only applied to non-NA below)

    } else if (na_action == "remove") {
      # Remove rows with missing factors
      rows_to_keep <- !DF_converted[[taxon_name]] %in% missing_taxa
      DF_converted <- DF_converted[rows_to_keep, ]
      DF_merged <- DF_merged[rows_to_keep, ]
      message("Removed ", length(missing_taxa), " taxa without conversion factors")

    } else if (na_action == "mean") {
      mean_factor <- mean(DF_merged[[factor_col]], na.rm = TRUE)
      DF_merged[[factor_col]][is.na(DF_merged[[factor_col]])] <- mean_factor
      message("Using mean factor (", round(mean_factor, 4),
              ") for taxa without specific factors")
    }
  }

  # Apply conversion ONLY to rows with valid factors
  for (col in biomass_cols) {
    # Identify rows with valid conversion factors
    has_factor <- !is.na(DF_merged[[factor_col]])

    # Apply conversion only to those rows
    DF_converted[has_factor, col] <- DF_converted[has_factor, col] *
      DF_merged[[factor_col]][has_factor]
  }

  return(DF_converted)
}


#' Get Internal Conversion Factors Database
#'
#' Returns the internal database of biomass conversion factors for benthic taxa.
#'
#' @return Data frame with conversion factors
#'
#' @details
#' This is an internal function that stores the conversion factors database.
#' Users typically don't need to call this directly - it's used automatically
#' by ConvertBiomass() when factors = NULL.
#'
#' @keywords internal
get_conversion_factors <- function() {

  # Build conversion factors database from literature
  factors <- data.frame(
    Group = c(
      # Echinodermata
      "Asteroidea", "Echinoidea", "Irregularia", "Ophiuroidea",
      "Holothuroidea", "Dendrochirota", "Aspidochirota", "Elasipoda", "Crinoidea",

      # Crustacea
      "Amphipoda", "Isopoda", "Decapoda", "Natantia", "Reptantia",
      "Mysidacea", "Cumacea", "Euphausiacea", "Cirripedia",

      # Mollusca
      "Bivalvia", "Pectinacea", "Tellinacea", "Mytilus edulis",
      "Streptoneura", "Ophistobranchia", "Nudibranchia", "Polyplacophora",
      "Benth. Cephalopoda", "Pelag. Cephalopoda",

      # Annelida
      "Oligochaeta", "Polych. Errantia", "Polych. Sedentaria",

      # Misc Groups
      "Porifera", "Actinaria", "Gorgonacea", "Anthozoa", "Bryozoa",
      "Nemertea", "Priapulida", "Sipunculida", "Ascidiae",

      # Insecta Larvae
      "Chironomidae", "Ephemeroptera", "Odonata", "Trichoptera",

      # Fish
      "Demersal Fish", "Pelagic Fish", "Mesopelagic Fish"
    ),

    WM_to_DM = c(
      # Echinodermata
      0.283, 0.333, 0.194, 0.46, 0.11, NA, 0.081, 0.108, 0.432,

      # Crustacea
      0.2, 0.2, 0.258, 0.267, 0.258, 0.197, 0.173, 0.254, 0.066,

      # Mollusca (SFDM)
      0.087, NA, NA, NA, 0.099, 0.077, 0.25, NA, 0.203, 0.203,

      # Annelida
      0.174, 0.199, 0.188,

      # Misc Groups
      0.186, 0.16, NA, NA, 0.199, 0.208, 0.095, 0.177, 0.063,

      # Insecta Larvae
      NA, NA, 0.226, NA,

      # Fish
      NA, NA, 0.162
    ),

    DM_to_AFDM = c(
      # Echinodermata
      0.438, 0.165, 0.121, 0.211, 0.476, 0.565, 0.499, 0.459, 0.238,

      # Crustacea
      0.72, 0.64, 0.68, 0.876, 0.68, 0.824, 0.63, 0.883, 0.79,

      # Mollusca (SFDM)
      0.831, 0.845, 0.833, NA, 0.838, 0.766, 0.693, NA, 0.9, 0.9,

      # Annelida
      0.323, 0.813, 0.732,

      # Misc Groups
      0.372, 0.864, 0.16, NA, 0.402, 0.816, 0.861, 0.654, 0.358,

      # Insecta Larvae
      0.931, 0.847, 0.888, 0.892,

      # Fish
      NA, NA, NA
    ),

    WM_to_AFDM = c(
      # Echinodermata
      0.124, 0.049, 0.023, 0.09, 0.112, 0.112, 0.04, 0.05, 0.08,

      # Crustacea
      0.16, 0.142, 0.18, 0.234, 0.18, 0.155, 0.075, 0.224, 0.039,

      # Mollusca (SFDM)
      0.057, 0.15, NA, NA, 0.076, 0.137, 0.173, 0.272, 0.2, 0.226,

      # Annelida
      NA, 0.169, 0.145,

      # Misc Groups
      0.075, 0.138, NA, NA, 0.08, 0.211, 0.065, 0.111, 0.023,

      # Insecta Larvae
      NA, NA, NA, NA,

      # Fish
      0.251, 0.207, 0.132
    ),

    J_per_mgDM = c(
      # Echinodermata
      9.11, 3.4, 1.66, 4.6, 11.27, 11.27, 11.27, 11.27, 5.1,

      # Crustacea
      NA, NA, NA, 16.23, NA, NA, NA, NA, NA,

      # Mollusca (SFDM)
      NA, 20.22, 18.47, 21.15, NA, NA, NA, NA, 20.4, 20.4,

      # Annelida
      NA, NA, NA,

      # Misc Groups
      7.75, NA, NA, NA, NA, 22.9, NA, NA, NA,

      # Insecta Larvae
      21.83, 22.07, 20.99, 21.52,

      # Fish
      NA, NA, NA
    ),

    J_per_mgAFDM = c(
      # Echinodermata
      20.81, 20.53, 20.53, 21.75, 22.95, 22.95, 22.95, 22.95, 21.44,

      # Crustacea
      22.74, 22.74, 22.26, 22.25, 22.26, 23, 22.74, 22.74, 22.74,

      # Mollusca (SFDM)
      22.79, 21.96, 22.18, 23.26, 23.63, 23.99, 23.27, 23.27, 22.03, 23.34,

      # Annelida
      23.33, 23.33, 23.33,

      # Misc Groups
      24.99, 21.54, 27.37, 24.46, 23.09, 23.33, 23.33, 23.33, 19.01,

      # Insecta Larvae
      23.44, 26.07, 23.65, 24.12,

      # Fish
      25.57, 23.32, 23.32
    ),

    stringsAsFactors = FALSE
  )

  return(factors)
}

#' View Available Conversion Factors
#'
#' Displays the complete internal database of biomass conversion factors,
#' showing which taxonomic groups have factors available for each conversion type.
#'
#' @param conversion Character. Optionally filter by conversion type:
#'   "WM_to_DM", "DM_to_AFDM", "WM_to_AFDM", "DM_to_Energy", or "AFDM_to_Energy".
#'   If NULL, shows all conversions (default = NULL).
#' @param group_pattern Character. Optional pattern to filter groups (e.g., "Amphi"
#'   to see all Amphipoda-related entries). Uses grep pattern matching.
#'
#' @return Data frame with available conversion factors.
#'
#' @details
#' This function helps identify the exact group names needed for ConvertBiomass().
#' Use it to check spelling and see which conversions are available for your taxa.
#'
#' @examples
#' \dontrun{
#' # View all available factors
#' view_conversion_factors()
#'
#' # View only WM to DM conversions
#' view_conversion_factors(conversion = "WM_to_DM")
#'
#' # Search for specific group
#' view_conversion_factors(group_pattern = "Amphi")
#'
#' # Search for all Crustacea
#' view_conversion_factors(group_pattern = "acea")
#'
#' # Check what's available for Mollusca
#' view_conversion_factors(group_pattern = "alv")
#' }
#'
#' @seealso \code{\link{ConvertBiomass}} for converting biomass values
#'
#' @export
view_conversion_factors <- function(conversion = NULL, group_pattern = NULL) {

  # Get internal database
  factors <- get_conversion_factors()

  # Filter by group pattern if provided
  if (!is.null(group_pattern)) {
    matching <- grep(group_pattern, factors$Group, ignore.case = TRUE)
    if (length(matching) == 0) {
      message("No groups matching pattern '", group_pattern, "' found.")
      return(invisible(NULL))
    }
    factors <- factors[matching, ]
  }

  # Filter by conversion type if provided
  if (!is.null(conversion)) {
    valid_conversions <- c("WM_to_DM", "DM_to_AFDM", "WM_to_AFDM",
                           "DM_to_Energy", "AFDM_to_Energy")
    if (!conversion %in% valid_conversions) {
      stop("Invalid conversion type. Must be one of: ",
           paste(valid_conversions, collapse = ", "))
    }

    # Map to column name
    factor_col <- switch(conversion,
                         "WM_to_DM" = "WM_to_DM",
                         "DM_to_AFDM" = "DM_to_AFDM",
                         "WM_to_AFDM" = "WM_to_AFDM",
                         "DM_to_Energy" = "J_per_mgDM",
                         "AFDM_to_Energy" = "J_per_mgAFDM"
    )

    # Show only groups with this conversion available
    factors_filtered <- factors[!is.na(factors[[factor_col]]), c("Group", factor_col)]

    cat("\nAvailable factors for", conversion, ":\n\n")
    print(factors_filtered, row.names = FALSE)
    cat("\nTotal:", nrow(factors_filtered), "groups\n")

    return(invisible(factors_filtered))
  }

  # Show all factors
  cat("\nComplete conversion factors database:\n\n")
  print(factors, row.names = FALSE)
  cat("\nTotal:", nrow(factors), "groups\n")
  cat("\nNote: NA values indicate conversion factor not available for that group.\n")
  cat("Use view_conversion_factors(group_pattern = 'search_term') to search.\n")

  return(invisible(factors))
}
