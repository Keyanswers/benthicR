#' Convert Biomass Units Using Taxon-Specific Conversion Factors
#'
#' @description
#' Converts biomass measurements between different units (e.g., wet mass to dry
#' mass, AFDM to energy) using taxon-specific conversion factors compiled from
#' extensive literature sources. The function includes a comprehensive internal
#' database of conversion factors for marine, and freshwater benthic taxa.
#'
#' @param DF Data frame with biomass data. Must contain:
#'   \itemize{
#'     \item A taxonomic column (matching taxon_col)
#'     \item One or more biomass columns to convert
#'   }
#' @param taxon_col Character or numeric. Column name or number containing
#'   taxonomic group (e.g., "Group", "Phylum", "Class"). Must match column
#'   names in the conversion factors table (default = "Group").
#' @param biomass_cols Numeric vector or character vector. Columns containing
#'   biomass values to convert. Can be column numbers (e.g., 3:5) or names
#'   (e.g., c("St1", "St2")).
#' @param conversion Character. Type of conversion to apply. Options:
#'   \itemize{
#'     \item \code{"WM_to_DM"} - Wet mass to dry mass
#'     \item \code{"DM_to_AFDM"} - Dry mass to ash-free dry mass
#'     \item \code{"WM_to_AFDM"} - Wet mass to ash-free dry mass (direct)
#'     \item \code{"DM_to_Energy"} - Dry mass to energy (J/mg DM)
#'     \item \code{"AFDM_to_Energy"} - AFDM to energy (J/mg AFDM)
#'     \item \code{"Water_percent"} - Water content as percentage
#'     \item \code{"Ash_percent"} - Ash content as percentage
#'     \item \code{"AFDW_percent"} - Ash-free dry weight as percentage
#'     \item \code{"Carbon_percent"} - Organic carbon as percentage
#'     \item \code{"WM_to_KCal"} - Wet mass to KCal (inland waters)
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
#'   \item \strong{DM_to_AFDM:} Dry mass → Ash-free dry mass (after combustion at 500°C)
#'   \item \strong{WM_to_AFDM:} Wet mass → AFDM (combined conversion)
#'   \item \strong{DM_to_Energy:} Dry mass → Energy content (Joules per mg DM)
#'   \item \strong{AFDM_to_Energy:} AFDM → Energy content (Joules per mg AFDM)
#'   \item \strong{Water_percent:} Water content as percentage of wet weight
#'   \item \strong{Ash_percent:} Ash content as percentage (of WW or DW depending on source)
#'   \item \strong{AFDW_percent:} Ash-free dry weight as percentage of wet weight
#'   \item \strong{Carbon_percent:} Organic carbon as percentage of wet weight
#'   \item \strong{WM_to_KCal:} Wet mass to KCal (for inland waters: fish = 1.0,
#'         zoobenthos/zooplankton = 0.83)
#' }
#'
#' \strong{Internal Database:}
#'
#' The function includes conversion factors for major benthic taxa compiled from
#' literature sources including:
#'
#' \strong{Marine Taxa:}
#' \itemize{
#'   \item \strong{Porifera:} Sponges
#'   \item \strong{Cnidaria:} Hydrozoa, Octocorallia, Scleractinia, Actiniaria, Anthozoa
#'   \item \strong{Nemertea:} Ribbon worms
#'   \item \strong{Annelida:} Polychaeta (Errantia, Sedentaria), Oligochaeta
#'   \item \strong{Sipuncula:} Peanut worms
#'   \item \strong{Mollusca:} Gastropoda, Bivalvia, Cephalopoda, Polyplacophora
#'   \item \strong{Arthropoda:} Amphipoda, Isopoda, Decapoda, Mysidacea, Cumacea,
#'         Euphausiacea, Cirripedia
#'   \item \strong{Echinodermata:} Asteroidea, Echinoidea, Ophiuroidea, Holothuroidea,
#'         Crinoidea
#'   \item \strong{Bryozoa:} Moss animals
#'   \item \strong{Tunicata:} Ascidians
#'   \item \strong{Pisces:} Demersal, Pelagic, Mesopelagic fishes
#' }
#'
#' \strong{Freshwater/Inland Taxa:}
#' \itemize{
#'   \item \strong{Insecta Larvae:} Chironomidae, Ephemeroptera, Odonata, Trichoptera
#'   \item \strong{Fish (inland):} Wet mass to dry mass (0.20), WM to KCal (1.0)
#'   \item \strong{Zoobenthos/Zooplankton (inland):} WM to DM (0.167), WM to KCal (0.83)
#' }
#'
#' \strong{Principal Literature Sources:}
#' \itemize{
#'   \item Vinogradov (1953) - Elementary chemical composition of marine organisms
#'   \item Waters (1977) - Secondary production in inland waters
#'   \item Brey et al. (1988) - Energy content of macrobenthic invertebrates
#'   \item Brey (2010, 2012) - Aquatic invertebrate respiration and production
#'   \item Ricciardi & Bourget (1998) - Weight-to-weight conversion factors
#'   \item Eleftheriou & Basford (1989) - Macrobenthic infauna of North Sea
#'   \item Rumohr et al. (1987) - Conversion factors for benthos
#'   \item Lie (1968) - Production studies in Puget Sound
#'   \item Cummins & Wuycheck (1971) - Caloric equivalents
#'   \item Benke et al. (1999) - Length-mass relationships for freshwater invertebrates
#' }
#'
#' \strong{Shell-Free Dry Mass (SFDM):}
#'
#' For molluscs, factors are based on shell-free dry mass. If your data includes
#' shells, pre-process to remove shell weight or adjust factors accordingly.
#'
#' @note
#' \itemize{
#'   \item Conversion factors are approximations and can vary with season, location,
#'         and individual size. For critical applications, consider using taxon-specific
#'         or site-specific factors.
#'   \item For inland waters, use Waters (1977) factors: WM_to_DM = 0.20 (fish) or
#'         0.167 (invertebrates).
#' }
#'
#' @seealso
#' \code{\link{view_conversion_factors}} to explore available conversion factors.
#'
#' @examples
#' \dontrun{
#' # Example 1: Convert wet mass to dry mass (marine)
#' biomass_data <- data.frame(
#'   Group = c("Amphipoda", "Bivalvia", "Polychaeta", "Asteroidea"),
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
#'   biomass_cols = 3:5,
#'   conversion = "WM_to_DM"
#' )
#'
#' print(biomass_DM)
#'
#' # Example 2: Convert dry mass to energy content
#' colnames(biomass_DM)[3:5] <- c("St1_DM", "St2_DM", "St3_DM")
#'
#' biomass_energy <- ConvertBiomass(
#'   DF = biomass_DM,
#'   taxon_col = "Group",
#'   biomass_cols = 3:5,
#'   conversion = "DM_to_Energy"
#' )
#'
#' # Example 3: Inland waters conversion (Waters 1977)
#' inland_data <- data.frame(
#'   Group = c("Fish Inland", "Zoobenthos Inland", "Zooplankton Inland"),
#'   Site1 = c(50.0, 25.5, 10.2),
#'   Site2 = c(45.3, 30.1, 12.5)
#' )
#'
#' inland_DM <- ConvertBiomass(
#'   DF = inland_data,
#'   taxon_col = "Group",
#'   biomass_cols = 2:3,
#'   conversion = "WM_to_DM"
#' )
#'
#' inland_KCal <- ConvertBiomass(
#'   DF = inland_data,
#'   taxon_col = "Group",
#'   biomass_cols = 2:3,
#'   conversion = "WM_to_KCal"
#' )
#'
#' # Example 4: Use custom conversion factors
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
#' # Example 5: Search for available groups
#' view_conversion_factors(group_pattern = "acea")      # Crustacea
#' view_conversion_factors(group_pattern = "Polych")    # Polychaeta
#' view_conversion_factors(group_pattern = "Echino")    # Echinodermata
#' view_conversion_factors(conversion = "Water_percent") # Water content
#'
#' # Example 6: Handle missing taxa
#' biomass_with_na <- data.frame(
#'   Group = c("Amphipoda", "UnknownGroup", "Bivalvia"),
#'   Species = c("Ampelisca", "Mystery", "Nucula"),
#'   St1 = c(10, 15, 20)
#' )
#'
#' result <- ConvertBiomass(
#'   DF = biomass_with_na,
#'   taxon_col = "Group",
#'   biomass_cols = 3,
#'   conversion = "WM_to_DM",
#'   na_action = "mean"
#' )
#' }
#' @export
#' @references
#' Barnes, A.T., Quetin, L.B. & Childress, J.J. (1976). Deep-sea macroplanktonic
#' sea cucumbers: suspended sediment feeders captured from deep submergence
#' vehicle. Science, 194, 1083-1085.
#'
#' Billett, D.S.M. (1991). Deep-sea holothurians. Oceanography and Marine
#' Biology: An Annual Review, 29, 259-317.
#'
#' Brey, T. (2010). An empirical model for estimating aquatic invertebrate
#' respiration. Methods in Ecology and Evolution, 1, 92-101.
#' doi: 10.1111/j.2041-210X.2009.00008.x
#'
#' Brey, T. (2012). A multi-parameter artificial neural network model to
#' estimate macrobenthic invertebrate productivity and production. Limnology
#' and Oceanography: Methods, 10, 581-589. doi: 10.4319/lom.2012.10.581
#'
#' Brey, T., Rumohr, H. & Ankar, S. (1988). Energy content of macrobenthic
#' invertebrates: general conversion factors from weight to energy. Journal
#' of Experimental Marine Biology and Ecology, 117, 271-278.
#' doi: 10.1016/0022-0981(88)90062-7
#'
#' Eleftheriou, A. & Basford, D.J. (1989). The macrobenthic infauna of the
#' offshore northern North Sea. Journal of the Marine Biological Association
#' of the United Kingdom, 69, 123-143. doi: 10.1017/S0025315400049158
#'
#' Galeron, J., Sibuet, M., Mahaut, M.-L. & Dinet, A. (2000). Variation in
#' structure and biomass of the benthic communities at three contrasting
#' sites in the tropical Northeast Atlantic. Marine Ecology Progress Series,
#' 197, 121-137.
#'
#' Ricciardi, A. & Bourget, E. (1998). Weight-to-weight conversion factors
#' for marine benthic macroinvertebrates. Marine Ecology Progress Series,
#' 163, 245-251.
#'
#' Smith, C.R. & Hamilton, S.C. (1983). Epibenthic megafauna of a bathyal
#' basin off southern California: patterns of abundance, biomass, and
#' dispersion. Deep-Sea Research, 30, 907-928.
#'
#' Steimle, F.W. & Terranova, R.T. (1985). Energetic equivalents of marine
#' organisms from the continental shelf of the temperate Northwest Atlantic.
#' Journal of Northwest Atlantic Fisheries Science, 6, 117-124.
#'
#' Torres, J.J., Belman, B.W. & Childress, J.J. (1979). Oxygen consumption
#' rates of midwater fishes as a function of depth of occurrence. Deep-Sea
#' Research, 26A, 185-197.
#'
#' Vinogradov, A.P. (1953). The elementary chemical composition of marine
#' organisms. Sears Foundation for Marine Research, Yale University, New Haven.
#'
#' Wacasey, J.W. & Atkinson, E.G. (1987). Energy values of marine benthic
#' invertebrates from the Canadian Arctic. Marine Ecology Progress Series,
#' 39, 243-250.
#'
#' Walker, M., Tyler, P.A. & Billett, D.S.M. (1987). Organic and calorific
#' content of body tissues of deep-sea elasipodid holothurians in the
#' northeast Atlantic Ocean. Marine Biology, 96, 277-282.
#'
#' Zibrowius, H. (1985). Scleractiniaires bathyaux et abyssaux de
#' l'Atlantique nord-oriental: campagnes Biogas (Polygas) et Incal. In:
#' Laubier, L. & Monniot, C. (eds.) Peuplements profonds du Golfe de
#' Gascogne. IFREMER, Paris, pp. 311-324.
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

  # Input validation (unchanged)
  if (!is.data.frame(DF)) {
    stop("DF must be a data frame")
  }

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

  if (is.character(biomass_cols)) {
    biomass_cols_nums <- which(colnames(DF) %in% biomass_cols)
    if (length(biomass_cols_nums) == 0) {
      stop("No matching biomass columns found in DF")
    }
    biomass_cols <- biomass_cols_nums
  }

  valid_conversions <- c("WM_to_DM", "DM_to_AFDM", "WM_to_AFDM",
                         "DM_to_Energy", "AFDM_to_Energy",
                         "Water_percent", "Ash_percent",
                         "AFDW_percent", "Carbon_percent", "WM_to_KCal")

  if (!conversion %in% valid_conversions) {
    stop("Invalid conversion type. Must be one of: ",
         paste(valid_conversions, collapse = ", "))
  }

  if (!na_action %in% c("warning", "remove", "mean")) {
    stop("na_action must be 'warning', 'remove', or 'mean'")
  }

  # Load factors
  if (is.null(factors)) {
    factors <- get_conversion_factors()
  } else {
    if (!is.data.frame(factors)) {
      stop("factors must be a data frame")
    }
  }

  # Validate factor column exists
  if (!conversion %in% colnames(factors)) {
    stop("Conversion factor '", conversion, "' not found in factors table")
  }

  # Merge and convert
  DF_converted <- DF
  DF$row_index <- seq_len(nrow(DF))

  DF_merged <- merge(DF, factors[, c("Group", conversion)],
                     by.x = taxon_name,
                     by.y = "Group",
                     all.x = TRUE)

  DF_merged <- DF_merged[order(DF_merged$row_index), ]
  DF_merged$row_index <- NULL
  DF_converted$row_index <- NULL

  # Handle missing factors
  missing_taxa <- unique(DF[[taxon_name]][is.na(DF_merged[[conversion]])])

  if (length(missing_taxa) > 0) {
    if (na_action == "warning") {
      warning("No conversion factors found for: ",
              paste(missing_taxa, collapse = ", "),
              ". Original values retained.")
    } else if (na_action == "remove") {
      rows_to_keep <- !DF_converted[[taxon_name]] %in% missing_taxa
      DF_converted <- DF_converted[rows_to_keep, ]
      DF_merged <- DF_merged[rows_to_keep, ]
      message("Removed ", length(missing_taxa), " taxa without conversion factors")
    } else if (na_action == "mean") {
      mean_factor <- mean(DF_merged[[conversion]], na.rm = TRUE)
      DF_merged[[conversion]][is.na(DF_merged[[conversion]])] <- mean_factor
      message("Using mean factor (", round(mean_factor, 4),
              ") for taxa without specific factors")
    }
  }

  # Apply conversion
  for (col in biomass_cols) {
    has_factor <- !is.na(DF_merged[[conversion]])
    DF_converted[has_factor, col] <- DF_converted[has_factor, col] *
      DF_merged[[conversion]][has_factor]
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

  # Build comprehensive conversion factors database from literature
  # Total: 87 taxonomic groups

  factors <- data.frame(
    Group = c(
      # ==================== PORIFERA ====================
      "Porifera",

      # ==================== CNIDARIA ====================
      "Hydrozoa", "Octocorallia", "Gorgonacea", "Zoantharia",
      "Actiniaria", "Scleractinia", "Anthozoa",

      # ==================== NEMERTEA ====================
      "Nemertea",

      # ==================== ANNELIDA ====================
      "Polychaeta", "Polychaeta Errantia", "Polychaeta Sedentaria",
      "Oligochaeta", "Annelida",

      # ==================== SIPUNCULA ====================
      "Sipuncula", "Sipunculida",

      # ==================== MOLLUSCA ====================
      "Aplacophora", "Gastropoda", "Gastropoda Streptoneura",
      "Gastropoda Opisthobranchia", "Gastropoda Nudibranchia",
      "Scaphopoda", "Bivalvia", "Bivalvia Pectinacea", "Bivalvia Tellinacea",
      "Mytilus edulis", "Cephalopoda", "Cephalopoda Benthica",
      "Cephalopoda Pelagica", "Polyplacophora", "Mollusca",

      # ==================== ARTHROPODA ====================
      "Pycnogonida", "Cirripedia", "Cumacea", "Tanaidacea", "Isopoda",
      "Amphipoda", "Mysidacea", "Euphausiacea", "Decapoda",
      "Decapoda Natantia", "Decapoda Reptantia", "Crustacea",

      # ==================== ECHINODERMATA ====================
      "Holothuroidea", "Holothuroidea translucida", "Holothuroidea densa",
      "Holothuroidea Dendrochirota", "Holothuroidea Aspidochirota",
      "Holothuroidea Elasipoda", "Asteroidea", "Ophiuroidea",
      "Echinoidea", "Echinoidea Regularia", "Echinoidea Irregularia",
      "Crinoidea", "Echinodermata",

      # ==================== BRYOZOA ====================
      "Bryozoa",

      # ==================== PRIAPULIDA ====================
      "Priapulida",

      # ==================== TUNICATA ====================
      "Tunicata", "Ascidiae",

      # ==================== INSECTA LARVAE ====================
      "Chironomidae", "Ephemeroptera", "Odonata", "Plecoptera",
      "Trichoptera", "Coleoptera", "Diptera", "Insecta Larvae",

      # ==================== PISCES ====================
      "Fish Demersal", "Fish Pelagic", "Fish Mesopelagic",
      "Fish Inland", "Zoobenthos Inland", "Zooplankton Inland"
    ),

    # ========== WATER CONTENT (% of WW) ==========
    Water_percent = c(
      # Porifera (1)
      80,
      # Cnidaria (7)
      85, 80, NA, NA, 85.6, NA, NA,
      # Nemertea (1)
      79.3,
      # Annelida (5)
      78.1, NA, NA, NA, NA,
      # Sipuncula (2)
      82, NA,
      # Mollusca (15)
      85, 75.4, NA, NA, NA, 85, 48, NA, NA, NA, 77, NA, NA, NA, NA,
      # Arthropoda (12)
      75, 30.4, 80.3, 80.3, 80.3, 78, NA, NA, NA, 76, 80.5, NA,
      # Echinodermata (13)
      NA, 96.2, 79, NA, NA, NA, 65, 63, 73.3, NA, 55, NA, NA,
      # Bryozoa (1)
      NA,
      # Priapulida (1)
      NA,
      # Tunicata (2)
      90, NA,
      # Insecta (8)
      NA, NA, NA, NA, NA, NA, NA, NA,
      # Pisces (6)
      NA, NA, NA, NA, NA, NA
    ),

    # ========== ASH CONTENT (% of DW unless noted) ==========
    Ash_percent = c(
      # Porifera (1)
      10,
      # Cnidaria (7)
      1.3, 10, NA, NA, 2.2, NA, NA,
      # Nemertea (1)
      2.6,
      # Annelida (5)
      4.4, 15.5, NA, NA, NA,
      # Sipuncula (2)
      7.8, NA,
      # Mollusca (15)
      2.5, 3, NA, NA, NA, 2.5, 70, NA, NA, NA, 1.5, NA, NA, 8.5, NA,
      # Arthropoda (12)
      5, 64, 27.9, 27.9, 27.9, NA, NA, NA, NA, 6.3, 3.6, 22.5,
      # Echinodermata (13)
      NA, 7.3, 3, NA, NA, NA, 50, 15, 25, NA, 85.3, NA, 8,
      # Bryozoa (1)
      NA,
      # Priapulida (1)
      NA,
      # Tunicata (2)
      93.1, NA,
      # Insecta (8)
      NA, NA, NA, NA, NA, NA, NA, NA,
      # Pisces (6)
      NA, NA, NA, NA, NA, NA
    ),

    # ========== ASH-FREE DRY WEIGHT (% of WW) ==========
    AFDW_percent = c(
      # Porifera (1)
      10,
      # Cnidaria (7)
      13.7, 10, NA, NA, 12.2, NA, NA,
      # Nemertea (1)
      18.1,
      # Annelida (5)
      17.5, NA, NA, NA, NA,
      # Sipuncula (2)
      10.2, NA,
      # Mollusca (15)
      12.5, 5, NA, NA, NA, 12.5, 15.4, NA, NA, NA, 21.5, NA, NA, NA, NA,
      # Arthropoda (12)
      20, 5.6, 14.5, 16.4, 18.3, 15.7, NA, NA, NA, 20.6, 12.2, NA,
      # Echinodermata (13)
      NA, 0.8, 10.4, NA, NA, NA, 20, 12, 3.9, NA, 1.3, NA, NA,
      # Bryozoa (1)
      NA,
      # Priapulida (1)
      NA,
      # Tunicata (2)
      3.5, NA,
      # Insecta (8)
      NA, NA, NA, NA, NA, NA, NA, NA,
      # Pisces (6)
      NA, NA, NA, NA, NA, NA
    ),

    # ========== ORGANIC CARBON (% of WW) ==========
    Carbon_percent = c(
      # Porifera (1)
      5.2,
      # Cnidaria (7)
      7.1, 5.2, NA, 2.1, 6.3, NA, NA,
      # Nemertea (1)
      9.3,
      # Annelida (5)
      9, NA, NA, NA, NA,
      # Sipuncula (2)
      5.3, NA,
      # Mollusca (15)
      6.5, 2.6, NA, NA, NA, 6.5, 8, NA, NA, NA, 11.1, NA, NA, NA, NA,
      # Arthropoda (12)
      10.4, 2.9, 7.5, 8.5, 9.5, 8.1, NA, NA, NA, 10.6, 6.4, NA,
      # Echinodermata (13)
      NA, 0.5, 5.4, NA, NA, NA, 10.4, 6.4, 2, NA, 0.7, NA, NA,
      # Bryozoa (1)
      NA,
      # Priapulida (1)
      NA,
      # Tunicata (2)
      1.7, NA,
      # Insecta (8)
      NA, NA, NA, NA, NA, NA, NA, NA,
      # Pisces (6)
      NA, NA, NA, NA, NA, NA
    ),

    # ========== WM to DM CONVERSION ==========
    WM_to_DM = c(
      # Porifera (1)
      0.186,
      # Cnidaria (7)
      NA, NA, NA, NA, 0.16, NA, NA,
      # Nemertea (1)
      0.208,
      # Annelida (5)
      NA, 0.199, 0.188, 0.174, 0.187,
      # Sipuncula (2)
      0.177, 0.177,
      # Mollusca (15)
      NA, NA, 0.099, 0.077, 0.25, NA, 0.087, NA, NA, NA, 0.203, 0.203, 0.203, NA, NA,
      # Arthropoda (12)
      NA, 0.066, 0.173, NA, 0.2, 0.2, 0.197, 0.254, 0.258, 0.267, 0.258, 0.226,
      # Echinodermata (13)
      0.11, NA, NA, NA, 0.081, 0.108, 0.283, 0.46, 0.333, 0.333, 0.194, 0.432, NA,
      # Bryozoa (1)
      0.199,
      # Priapulida (1)
      0.095,
      # Tunicata (2)
      0.063, 0.063,
      # Insecta (8)
      NA, NA, 0.226, NA, NA, NA, NA, 0.21,
      # Pisces (6)
      NA, NA, 0.162, 0.20, 0.167, 0.167
    ),

    # ========== DM to AFDM CONVERSION ==========
    DM_to_AFDM = c(
      # Porifera (1)
      0.372,
      # Cnidaria (7)
      NA, NA, 0.16, NA, 0.864, NA, NA,
      # Nemertea (1)
      0.816,
      # Annelida (5)
      NA, 0.813, 0.732, 0.323, 0.623,
      # Sipuncula (2)
      0.654, 0.654,
      # Mollusca (15)
      NA, NA, 0.838, 0.766, 0.693, NA, 0.831, 0.845, 0.833, NA, 0.9, 0.9, 0.9, NA, 0.801,
      # Arthropoda (12)
      NA, 0.79, 0.63, NA, 0.64, 0.72, 0.824, 0.883, 0.68, 0.876, 0.68, 0.742,
      # Echinodermata (13)
      0.476, NA, NA, 0.565, 0.499, 0.459, 0.438, 0.211, 0.165, 0.165, 0.121, 0.238, NA,
      # Bryozoa (1)
      0.402,
      # Priapulida (1)
      0.861,
      # Tunicata (2)
      0.358, 0.358,
      # Insecta (8)
      0.931, 0.847, 0.888, NA, 0.892, NA, NA, 0.942,
      # Pisces (6)
      NA, NA, NA, NA, 0.9, 0.9
    ),

    # ========== WM to AFDM CONVERSION ==========
    WM_to_AFDM = c(
      # Porifera (1)
      0.075,
      # Cnidaria (7)
      NA, NA, NA, NA, 0.138, NA, NA,
      # Nemertea (1)
      0.211,
      # Annelida (5)
      0.132, 0.169, 0.145, NA, 0.157,
      # Sipuncula (2)
      0.111, 0.111,
      # Mollusca (15)
      NA, NA, 0.076, 0.137, 0.173, NA, 0.057, 0.15, NA, NA, 0.2, 0.226, 0.2, 0.272, 0.144,
      # Arthropoda (12)
      NA, 0.039, 0.075, NA, 0.142, 0.16, 0.155, 0.224, 0.18, 0.234, 0.18, 0.169,
      # Echinodermata (13)
      0.112, NA, NA, 0.112, 0.04, 0.05, 0.124, 0.09, 0.049, 0.049, 0.023, 0.08, NA,
      # Bryozoa (1)
      0.08,
      # Priapulida (1)
      0.065,
      # Tunicata (2)
      0.023, 0.023,
      # Insecta (8)
      NA, NA, NA, NA, NA, NA, NA, NA,
      # Pisces (6)
      0.251, 0.207, 0.132, NA, 1.11, 1.11
    ),

    # ========== ENERGY CONTENT (J/mg DM) ==========
    J_per_mgDM = c(
      # Porifera (1)
      7.75,
      # Cnidaria (7)
      NA, NA, NA, NA, NA, NA, NA,
      # Nemertea (1)
      22.9,
      # Annelida (5)
      NA, NA, NA, NA, NA,
      # Sipuncula (2)
      NA, NA,
      # Mollusca (15)
      NA, NA, NA, NA, NA, NA, NA, 20.22, 18.47, 21.15, NA, NA, 20.4, NA, NA,
      # Arthropoda (12)
      NA, NA, NA, NA, NA, NA, NA, NA, NA, 16.23, NA, NA,
      # Echinodermata (13)
      11.27, NA, NA, 11.27, 11.27, 11.27, 9.11, 4.6, 3.4, 3.4, 1.66, 5.1, NA,
      # Bryozoa (1)
      NA,
      # Priapulida (1)
      NA,
      # Tunicata (2)
      NA, NA,
      # Insecta (8)
      21.83, 22.07, 20.99, NA, 21.52, NA, NA, 22.44,
      # Pisces (6)
      NA, NA, NA, NA, NA, NA
    ),

    # ========== ENERGY CONTENT (J/mg AFDM) ==========
    J_per_mgAFDM = c(
      # Porifera (1)
      24.99,
      # Cnidaria (7)
      NA, NA, 27.37, NA, 21.54, NA, 24.46,
      # Nemertea (1)
      23.33,
      # Annelida (5)
      23.33, 23.33, 23.33, 23.33, 23.33,
      # Sipuncula (2)
      23.33, 23.33,
      # Mollusca (15)
      NA, NA, 23.63, 23.99, 23.27, NA, 22.79, 21.96, 22.18, 23.26, 22.03, 23.34, 22.03, 23.27, 23.04,
      # Arthropoda (12)
      NA, 22.74, 22.74, NA, 22.74, 22.74, 23.00, 22.74, 22.26, 22.25, 22.26, 22.57,
      # Echinodermata (13)
      22.95, NA, NA, 22.95, 22.95, 22.95, 20.81, 21.75, 20.53, 20.53, 20.53, 21.44, NA,
      # Bryozoa (1)
      23.09,
      # Priapulida (1)
      23.33,
      # Tunicata (2)
      19.01, 19.01,
      # Insecta (8)
      23.44, 26.07, 23.65, NA, 24.12, NA, NA, 23.81,
      # Pisces (6)
      25.57, 23.32, 23.32, 25.57, NA, NA
    ),

    # ========== WM to KCal CONVERSION ==========
    WM_to_KCal = c(
      # Porifera (1)
      NA,
      # Cnidaria (7)
      NA, NA, NA, NA, NA, NA, NA,
      # Nemertea (1)
      NA,
      # Annelida (5)
      NA, NA, NA, NA, NA,
      # Sipuncula (2)
      NA, NA,
      # Mollusca (15)
      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      # Arthropoda (12)
      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      # Echinodermata (13)
      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      # Bryozoa (1)
      NA,
      # Priapulida (1)
      NA,
      # Tunicata (2)
      NA, NA,
      # Insecta (8)
      NA, NA, NA, NA, NA, NA, NA, NA,
      # Pisces (6)
      NA, NA, NA, 1.0, 0.83, 0.83
    ),

    # ========== SOURCE REFERENCES ==========
    Source = c(
      # Porifera (1)
      "Vinogradov (1953); Brey et al. refs",
      # Cnidaria (7)
      "Vinogradov (1953)", "Vinogradov (1953)", "Brey et al. refs",
      "Galeron (unpubl.)", "Galeron (unpubl.); Brey et al. refs",
      "Zibrowius (1985)", "Brey et al. refs",
      # Nemertea (1)
      "Vinogradov (1953); Brey et al. refs",
      # Annelida (5)
      "Vinogradov (1953); Rumohr et al. (1987); Ricciardi & Bourget (1998); Lie (1968)",
      "Eleftheriou & Basford (1989); Brey et al. refs",
      "Brey et al. refs", "Brey et al. refs", "Brey et al. refs",
      # Sipuncula (2)
      "Vinogradov (1953)", "Brey et al. refs",
      # Mollusca (15)
      "Vinogradov (1953)", "Vinogradov (1953)", "Brey et al. refs",
      "Brey et al. refs", "Brey et al. refs", "Vinogradov (1953)",
      "Vinogradov (1953); Brey et al. refs", "Brey et al. refs",
      "Brey et al. refs", "Brey et al. refs", "Vinogradov (1953)",
      "Brey et al. refs", "Brey et al. refs", "Brey et al. refs",
      "Eleftheriou & Basford (1989); Brey et al. refs",
      # Arthropoda (12)
      "Vinogradov (1953)", "Vinogradov (1953); Brey et al. refs",
      "Vinogradov (1953); Brey et al. refs", "Vinogradov (1953)",
      "Vinogradov (1953); Brey et al. refs", "Vinogradov (1953); Brey et al. refs",
      "Brey et al. refs", "Brey et al. refs", "Brey et al. refs",
      "Vinogradov (1953); Brey et al. refs", "Vinogradov (1953); Brey et al. refs",
      "Eleftheriou & Basford (1989); Brey et al. refs",
      # Echinodermata (13)
      "Brey et al. refs", "Barnes et al. (1976)", "Vinogradov (1953)",
      "Brey et al. refs", "Brey et al. refs", "Brey et al. refs",
      "Vinogradov (1953); Brey et al. refs", "Smith & Hamilton (1983); Brey et al. refs",
      "Vinogradov (1953); Brey et al. refs", "Brey et al. refs",
      "Vinogradov (1953); Brey et al. refs", "Brey et al. refs",
      "Eleftheriou & Basford (1989)",
      # Bryozoa (1)
      "Brey et al. refs",
      # Priapulida (1)
      "Brey et al. refs",
      # Tunicata (2)
      "Vinogradov (1953)", "Brey et al. refs",
      # Insecta (8)
      "Brey et al. refs", "Brey et al. refs", "Brey et al. refs",
      "Waters (1977)", "Brey et al. refs", "Waters (1977)",
      "Waters (1977)", "Brey et al. refs",
      # Pisces (6)
      "Brey et al. refs", "Brey et al. refs", "Brey et al. refs",
      "Waters (1977)", "Waters (1977)", "Waters (1977)"
    ),

    stringsAsFactors = FALSE
  )

  return(factors)
}

#' View Available Conversion Factors
#'
#' Displays the internal database of biomass conversion factors compiled from
#' literature sources. Can filter by conversion type, taxonomic group pattern,
#' or both.
#'
#' @param conversion Character. Type of conversion to display. If NULL, shows
#'   all available conversions for matching groups. Options:
#'   \itemize{
#'     \item \code{"WM_to_DM"} - Wet mass to dry mass
#'     \item \code{"DM_to_AFDM"} - Dry mass to ash-free dry mass
#'     \item \code{"WM_to_AFDM"} - Wet mass to ash-free dry mass
#'     \item \code{"DM_to_Energy"} - Dry mass to energy (J/mg DM)
#'     \item \code{"AFDM_to_Energy"} - AFDM to energy (J/mg AFDM)
#'     \item \code{"Water_percent"} - Water content as percentage
#'     \item \code{"Ash_percent"} - Ash content as percentage
#'     \item \code{"AFDW_percent"} - Ash-free dry weight as percentage
#'     \item \code{"Carbon_percent"} - Organic carbon as percentage
#'     \item \code{"WM_to_KCal"} - Wet mass to KCal (inland waters)
#'   }
#' @param group_pattern Character. Regular expression pattern to filter
#'   taxonomic groups. Case-insensitive. If NULL, shows all groups (default = NULL).
#' @param show_sources Logical. If TRUE, includes literature source column in
#'   output (default = TRUE).
#'
#' @return Invisibly returns a data frame with the displayed conversion factors.
#'
#' @details
#' This function provides an interface to explore the internal conversion
#' factors database used by \code{\link{ConvertBiomass}}. The database includes
#' factors for over 80 taxonomic groups of marine and freshwater benthic
#' organisms.
#'
#' \strong{Search Examples:}
#' \itemize{
#'   \item \code{view_conversion_factors()} - Show all groups and conversions
#'   \item \code{view_conversion_factors(conversion = "WM_to_DM")} - Only WM→DM factors
#'   \item \code{view_conversion_factors(group_pattern = "Polych")} - Groups containing "Polych"
#'   \item \code{view_conversion_factors(conversion = "AFDM_to_Energy", group_pattern = "acea")} - AFDM energy for crustaceans
#' }
#'
#' \strong{Note:} Energy conversions use column names \code{J_per_mgDM} and
#' \code{J_per_mgAFDM} internally but can be requested as \code{"DM_to_Energy"}
#' and \code{"AFDM_to_Energy"} for clarity.
#'
#' @seealso
#' \code{\link{ConvertBiomass}} for applying conversion factors to data.
#'
#' @examples
#' \dontrun{
#' # View all available factors
#' view_conversion_factors()
#'
#' # View only wet mass to dry mass conversions
#' view_conversion_factors(conversion = "WM_to_DM")
#'
#' # Search for all Crustacea groups
#' view_conversion_factors(group_pattern = "acea")
#'
#' # Search for Polychaeta energy conversions
#' view_conversion_factors(conversion = "AFDM_to_Energy", group_pattern = "Polych")
#'
#' # Hide source column for cleaner output
#' view_conversion_factors(conversion = "WM_to_DM", show_sources = FALSE)
#' }
#' @export

view_conversion_factors <- function(conversion = NULL, group_pattern = NULL,
                                    show_sources = TRUE) {

  factors <- get_conversion_factors()

  # Filter by group pattern FIRST (if provided)
  if (!is.null(group_pattern)) {
    matching <- grep(group_pattern, factors$Group, ignore.case = TRUE)
    if (length(matching) == 0) {
      message("No groups matching pattern '", group_pattern, "' found.")
      return(invisible(NULL))
    }
    factors <- factors[matching, ]
  }

  # If NO conversion specified, show ALL available factors for filtered groups
  if (is.null(conversion)) {

    # Identify which conversion columns have at least one non-NA value
    conversion_cols <- c("WM_to_DM", "DM_to_AFDM", "WM_to_AFDM",
                         "J_per_mgDM", "J_per_mgAFDM",
                         "Water_percent", "Ash_percent",
                         "AFDW_percent", "Carbon_percent", "WM_to_KCal")

    # Keep only columns that have at least one non-NA value in the filtered data
    cols_with_data <- conversion_cols[sapply(factors[conversion_cols],
                                             function(x) any(!is.na(x)))]

    # Select columns to display: Group + available conversions + (optionally) Source
    if (length(cols_with_data) > 0) {
      cols_to_show <- c("Group", cols_with_data)
    } else {
      cols_to_show <- "Group"
    }

    if (show_sources) {
      cols_to_show <- c(cols_to_show, "Source")
    }

    cat("\n========================================\n")
    if (!is.null(group_pattern)) {
      cat("GROUPS MATCHING PATTERN: '", group_pattern, "'\n", sep = "")
    } else {
      cat("COMPLETE CONVERSION FACTORS DATABASE\n")
    }
    cat("========================================\n\n")
    cat("Total groups:", nrow(factors), "\n\n")

    print(factors[, cols_to_show, drop = FALSE], row.names = FALSE)

    cat("\n")
    if (!is.null(group_pattern)) {
      cat("Use view_conversion_factors(conversion = 'type', group_pattern = '",
          group_pattern, "') to see only specific conversion.\n", sep = "")
    } else {
      cat("Use view_conversion_factors(conversion = 'type') to see specific factors.\n")
      cat("Use view_conversion_factors(group_pattern = 'taxon') to search.\n")
    }

    return(invisible(factors))
  }

  # If conversion IS specified, validate and filter
  valid_conversions <- c("WM_to_DM", "DM_to_AFDM", "WM_to_AFDM",
                         "DM_to_Energy", "AFDM_to_Energy",
                         "Water_percent", "Ash_percent",
                         "AFDW_percent", "Carbon_percent", "WM_to_KCal")

  # Map friendly names to actual column names
  conversion_col <- switch(conversion,
                           "DM_to_Energy" = "J_per_mgDM",
                           "AFDM_to_Energy" = "J_per_mgAFDM",
                           conversion)

  if (!conversion %in% valid_conversions) {
    stop("Invalid conversion type. Must be one of: ",
         paste(valid_conversions, collapse = ", "))
  }

  # Show only groups with this conversion available
  factors_filtered <- factors[!is.na(factors[[conversion_col]]), ]

  if (nrow(factors_filtered) == 0) {
    message("No conversion factors available for: ", conversion)
    return(invisible(NULL))
  }

  # Select columns to display
  cols_to_show <- c("Group", conversion_col)
  if (show_sources) {
    cols_to_show <- c(cols_to_show, "Source")
  }

  cat("\n========================================\n")
  cat("Conversion factors for:", conversion)
  if (!is.null(group_pattern)) {
    cat(" (matching '", group_pattern, "')", sep = "")
  }
  cat("\n========================================\n\n")
  print(factors_filtered[, cols_to_show], row.names = FALSE)
  cat("\nTotal:", nrow(factors_filtered), "groups\n")

  return(invisible(factors_filtered))
}



