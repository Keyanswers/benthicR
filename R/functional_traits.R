#' Calculate Functional Trait Abundances from Fuzzy-Coded Data
#'
#' Aggregates species abundances by functional traits using fuzzy coding.
#' This function takes a fuzzy-coded trait matrix (where species can have
#' partial membership in multiple trait categories) and calculates total
#' abundance for each trait category.
#'
#' @param Tr Data frame with fuzzy-coded trait information. Must contain:
#'   \itemize{
#'     \item \code{Phylum}: Taxonomic phylum
#'     \item \code{Species}: Species name
#'     \item Trait columns: Numeric values 0-1 indicating trait affinity
#'   }
#' @param Dat Data frame with abundance data. Must contain:
#'   \itemize{
#'     \item \code{Phylum}: Taxonomic phylum (matching Tr)
#'     \item \code{Species}: Speciess name (matching Tr)
#'     \item Abundance columns: Numeric abundances per sample/station
#'   }
#' @param abundance_cols Numeric vector or character vector specifying which
#'   columns in Dat contain abundance data. Can be column numbers (e.g., 3:10)
#'   or column names (e.g., c("St1", "St2", "St3")).
#' @param threshold Numeric. Minimum trait affinity value (0-1) to include a
#'   species in that trait category. Speciess with affinity >= threshold are
#'   considered to possess the trait (default = 0.1).
#' @param trait_cols Character vector of trait column names to process. If NULL,
#'   processes all numeric columns in Tr except "Phylum" and "Species" (default = NULL).
#'
#' @return A data frame with:
#'   \itemize{
#'     \item Rows: Functional traits
#'     \item Columns: Samples/stations
#'     \item Values: Total abundance of species possessing each trait
#'   }
#'
#' @details
#' \strong{Fuzzy Coding Approach:}
#'
#' Fuzzy coding allows species to have partial membership in multiple trait
#' categories. For example, a species might be:
#' \itemize{
#'   \item 70\% predator (0.7)
#'   \item 30\% omnivore (0.3)
#'   \item 0\% herbivore (0.0)
#' }
#'
#' This function uses a **threshold approach**:
#' \itemize{
#'   \item Speciess with trait affinity >= threshold are included
#'   \item Their FULL abundance is counted for that trait
#'   \item A species can contribute to multiple traits if it exceeds threshold
#'     for each
#' }
#'
#' \strong{Processing Steps:}
#' \enumerate{
#'   \item For each trait, identify species with affinity >= threshold
#'   \item Merge these species with abundance data
#'   \item Sum abundances across all species possessing that trait
#'   \item Return matrix of trait abundances
#' }
#'
#' \strong{Common Functional Traits for Benthic Communities:}
#'
#' \emph{Feeding mode:}
#' \itemize{
#'   \item Predator, Herbivore, Omnivore
#'   \item Particles (suspension feeders)
#'   \item Sediment (deposit feeders)
#' }
#'
#' \emph{Mobility:}
#' \itemize{
#'   \item Motile, Sedentary, Sessile
#'   \item Hmotile (highly motile)
#' }
#'
#' \emph{Living habit:}
#' \itemize{
#'   \item Burrow, Surface, Attached
#'   \item Epifauna, Infauna
#' }
#'
#' \emph{Body size:}
#' \itemize{
#'   \item VoSmall, Small, Medium, Large, VLarge
#' }
#'
#' \emph{Longevity:}
#' \itemize{
#'   \item ALeast5 (<5 years)
#'   \item B5to10 (5-10 years)
#'   \item B10to50 (10-50 years)
#'   \item Older50 (>50 years)
#' }
#'
#' \emph{Fragility/Protection:}
#' \itemize{
#'   \item Hard, Strong, Fragile
#'   \item NoPrtec (no protection)
#' }
#'
#' @note
#' The threshold approach is simpler than weighted abundance (multiplying
#' abundance by trait affinity) but is commonly used in functional ecology.
#' If you need weighted abundances, consider using the FD or mFD packages.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with fuzzy-coded traits
#' # Trait matrix (fuzzy coding)
#' traits <- data.frame(
#'   Phylum = c("Arthropoda", "Arthropoda", "Arthropoda", "Echinodermata"),
#'   Species = c("Achelous floridanus", "Anisopagurus pygmaeus", "Hepatus pudibundus",
#'              "Astropecten alligator"),
#'   Predator = c(1, 0, 0, 1),
#'   Herbivore = c(0, 0, 0, 0),
#'   Omnivore = c(0, 0, 1, 0),
#'   PP = c(0, 0, 0, 0),
#'   Particles = c(0, 0, 0, 0),
#'   Sediment = c(0, 1, 0, 0),
#'   Hmotile = c(1, 0, 0, 0),
#'   Motile = c(0, 1, 0.75, 0.75),
#'   Sedentary = c(0, 0, 0, 0),
#'   Sessile = c(0, 0, 0, 0),
#'   Burrow = c(0, 0, 0, 0),
#'   Surface = c(1, 1, 0.25, 0.25),
#'   Attached = c(0, 0, 0.75, 0.75)
#' )
#'
#' print(traits)
#'
#'# Abundance data
#' abundance <- data.frame(
#'   Phylum = c("Arthropoda", "Arthropoda", "Arthropoda", "Echinodermata"),
#'   Species = c("Achelous floridanus", "Anisopagurus pygmaeus", "Hepatus pudibundus",
#'               "Astropecten alligator"),
#'   St1 = c(15, 25, 8, 5),
#'   St2 = c(12, 30, 6, 4),
#'   St3 = c(18, 22, 10, 6),
#'   St4 = c(10, 28, 5, 3),
#'   St5 = c(20, 35, 12, 7)
#' )
#'
#' print(abundance)
#'
#'# Calculate trait abundances using column numbers
#' trait_abund <- Functional(Tr = traits,
#'                            Dat = abundance,
#'                            abundance_cols = 3:7)
#'
#' trait_abund <- trait_abund[rowSums(trait_abund) > 0,]
#'
#' print(trait_abund)
#'
#'#            St1 St2 St3
#'# Predator   42  34  50
#'# Herbivore  20  25  18
#'# Omnivore    8   6  10
#'# Sediment   45  55  40
#'# Hmotile    15  12  18
#'# Motile     80  83  82
#'# Surface    95  95 100
#'# Attached   25  20  30
#'
#'# Example 2: Using column names
#' trait_abund <- Functional(Tr = traits,
#'                            Dat = abundance,
#'                            abundance_cols = c("St1", "St2", "St3", "St4","St5"),
#'                            threshold = 0.5)  # Stricter threshold
#'
#' trait_abund <-  trait_abund[rowSums(trait_abund) > 0,]
#' print(trait_abund)
#'
#'#            St1 St2 St3 St4 St5
#'# Predator   42  34  50  28  56
#'# Herbivore  20  25  18  22  28
#'# Omnivore    8   6  10   5  12
#'# Sediment   45  55  40  50  63
#'# Hmotile    15  12  18  10  20
#'# Motile     80  83  82  73 111
#'# Surface    70  75  70  67  97
#'# Attached   25  20  30  16  34
#'
#' # Example 3: Select specific traits only
#' trait_abund <- Functional(Tr = traits,
#'                           Dat = abundance,
#'                           abundance_cols = 3:5,
#'                           trait_cols = c("Predator", "Sediment", "Motile"))
#'
#' # Example 4: Visualize with RadChart
#' # First, need to transpose for RadChart (expects stations in rows)
#' trait_abund_t <- as.data.frame(t(trait_abund))
#'
#' trait_abund_t <- trait_abund_t[,colSums(trait_abund_t) > 0]
#' print(trait_abund_t)
#'
#'#        Predator Herbivore Omnivore Sediment Hmotile Motile Surface Attached
#'# St1       42        20        8       45      15     80      70       25
#'# St2       34        25        6       55      12     83      75       20
#'# St3       50        18       10       40      18     82      70       30
#'# St4       28        22        5       50      10     73      67       16
#'# St5       56        28       12       63      20    111      97       34
#'
#' RadChart(trait_abund_t, shape = "w", layout = c(2, 4))
#'
#' # Example 5: Real workflow
#' # Load your fuzzy-coded trait matrix
#' my_traits <- read.csv("fuzzy_traits.csv")
#'
#' # Load abundance data
#' my_abundance <- read.csv("abundance.csv")
#'
#' # Calculate trait abundances for stations 3-15
#' trait_results <- Functional(Tr = my_traits,
#'                             Dat = my_abundance,
#'                             abundance_cols = 3:15,
#'                             threshold = 0.1)
#'
#' # Visualize
#' trait_results_t <- as.data.frame(t(trait_results))
#' RadChart(trait_results_t, shape = "w")
#' }
#'
#' @references
#' Andrade, L. S., et al. (2016). The assemblage composition and structure of swimming crabs (Portunoidea)
#' in continental shelf waters of southeastern Brazil. Journal of the Marine Biological Association of the
#' United Kingdom, 96(5): 1021-1030.
#'
#' Biological Traits Information Catalogue (BIOTIC) - Marine Biological Association of the UK
#' https://www.marlin.ac.uk/biotic/
#'
#' Bitter, R. & Penchaszadeh, P. E. (1983). Ecología trófica de dos estrellas de mar del género Astropecten
#' coexistentes en Golfo Triste, Venezuela. Boletín del Instituto Oceanográfico de Venezuela, Universidad de Oriente, 22(1-2): 63-73.
#'
#' Bremner, J. (2008). "Speciess' traits and ecological functioning in marine conservation and management".
#' Journal of Experimental Marine Biology and Ecology, 366(1-2): 37-47.
#'
#' Chevene, F., Doléadec, S., & Chessel, D. (1994). A fuzzy coding approach
#' for the analysis of long-term ecological data. Freshwater Biology, 31(3), 295-309.
#'
#' Feder, H. (1963). Gastropod defensive responses and their effectiveness in reducing predation by
#' starfishes. Ecology, 44(3): 505-512.
#'
#' Hay, W. P. & Shore, C. A. (1905). The decapod crustaceans of Beaufort, N. C. and the surrounding region.
#' Bulletin of the Bureau of Fisheries, 23: 369-413. Available at: https://www.harteresearchinstitute.org/sites/default/files/inline-files/S72-S85.pdf (Accessed: 26 February 2019).
#'
#' Heck, K. L. (1977). Comparative species richness, composition, and abundance of invertebrates in Caribbean
#' seagrass (Thalassia testudinum) meadows (Panama). Marine Biology, 41: 335-348.
#'
#' Jesse, S. & Stotz, W. (2002). Spatio-temporal distribution patterns of the crab assemblage in the
#' shallow subtidal of the north Chilean Pacific coast. Crustaceana, 75(10): 1161-1200.
#'
#' Lemaitre, R. & McLaughlin, P. A. (1996). Revision of Pylopagurus and Tomopagurus (Crustacea: Decapoda:Paguridae),
#' with the descriptions of new genera and species. Part V. Anisopagurus McLaughlin, Manucomplanus McLaughlin, and
#' Protoniopagurus new genus. Bulletin of Marine Science, 59(1): 89-141.
#'
#' Medina-Mantelatto, F. L. & Petracco, M. (1997). Natural diet of the crab Hepatus pudibundus (Brachyura: Calappidae)
#' in Fortaleza Bay, Ubatuba (SP), Brazil. Journal of Crustacean Biology, 17(3): 440-446.
#'
#' O'Connor, W. A., et al. (2014). Difference in preference size of Perna viridis for optimal foraging of the sea star
#' Echinaster spinulosus. Journal of Shellfish Research, 33(2): 435-440.
#'
#' Telford, M., et al. (1987). Feeding activities of two species of Clypeaster (Echinoidea, Clypeasteroida):
#' Further evidence of clypeasteroid resource partitioning. Biological Bulletin, 172(3): 324-336.
#'
#' Willems, T. (2008). Impact of Guyana seabob trawl fishery on marine habitats and ecosystems:
#' A preliminary assessment. Report, 41 pp.
#'
#' WoRMS Editorial Board (2025). World Register of Marine Species. Available at: http://www.marinespecies.org (Accessed: 26 February 2019).
#'
#' @seealso \code{\link{RadChartMulti}} for visualizing trait abundances across groups,
#'   \code{\link{RadChart}} for comparing trait profiles,
#'   \code{\link{ExploreComm}} for exploring taxonomic composition
#'
#' @export

Functional <- function(Tr,
                       Dat,
                       abundance_cols,
                       threshold = 0.1,
                       trait_cols = NULL) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  if (!is.data.frame(Tr) || !is.data.frame(Dat)) {
    stop("Tr and Dat must be data frames")
  }

  if (!all(c("Phylum", "Species") %in% colnames(Tr))) {
    stop("Tr must contain columns 'Phylum' and 'Species'")
  }

  if (!all(c("Phylum", "Species") %in% colnames(Dat))) {
    stop("Dat must contain columns 'Phylum' and 'Species'")
  }

  if (missing(abundance_cols)) {
    stop("abundance_cols must be specified (column numbers or names)")
  }

  if (threshold < 0 || threshold > 1) {
    stop("threshold must be between 0 and 1")
  }

  # ============================================================================
  # PREPARE ABUNDANCE COLUMNS
  # ============================================================================

  # Convert abundance_cols to column numbers if they're names
  if (is.character(abundance_cols)) {
    abundance_cols_nums <- which(colnames(Dat) %in% abundance_cols)
    if (length(abundance_cols_nums) == 0) {
      stop("No matching abundance columns found in Dat")
    }
    abundance_cols <- abundance_cols_nums
  }

  # Validate abundance columns exist
  if (any(abundance_cols > ncol(Dat)) || any(abundance_cols < 1)) {
    stop("Invalid abundance column indices")
  }

  # ============================================================================
  # DETERMINE TRAIT COLUMNS
  # ============================================================================

  if (is.null(trait_cols)) {
    # All numeric columns except Phylum and Species
    trait_cols <- colnames(Tr)[!colnames(Tr) %in% c("Phylum", "Species")]
    trait_cols <- trait_cols[sapply(Tr[, trait_cols, drop = FALSE], is.numeric)]
  } else {
    # Validate provided trait columns exist
    missing_traits <- trait_cols[!trait_cols %in% colnames(Tr)]
    if (length(missing_traits) > 0) {
      stop("Trait columns not found in Tr: ", paste(missing_traits, collapse = ", "))
    }
  }

  if (length(trait_cols) == 0) {
    stop("No trait columns found in Tr")
  }

  # ============================================================================
  # INITIALIZE RESULTS MATRIX
  # ============================================================================

  sample_names <- colnames(Dat)[abundance_cols]
  results <- matrix(0, nrow = length(trait_cols), ncol = length(abundance_cols))
  rownames(results) <- trait_cols
  colnames(results) <- sample_names

  # ============================================================================
  # PROCESS EACH TRAIT
  # ============================================================================

  for (i in seq_along(trait_cols)) {
    trait <- trait_cols[i]

    # Filter species with this trait above threshold
    species_with_trait <- Tr[Tr[[trait]] >= threshold, c("Phylum", "Species")]

    if (nrow(species_with_trait) > 0) {
      # Merge with abundance data
      merged <- merge(Dat, species_with_trait, by = c("Phylum", "Species"))

      if (nrow(merged) > 0) {
        # Sum abundances for this trait
        results[i, ] <- colSums(merged[, abundance_cols, drop = FALSE], na.rm = TRUE)
      }
    }
  }

  # ============================================================================
  # CONVERT TO DATA FRAME AND RETURN
  # ============================================================================

  results_df <- as.data.frame(results)

  return(results_df)
}
