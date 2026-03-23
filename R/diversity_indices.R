#' Calculate Shannon Diversity Index
#'
#' Calculates the Shannon diversity index (Shannon-Wiener) for benthic community data
#' with flexible input formats and logarithm base options.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param base Logarithm base for calculation. Options are:
#'   \itemize{
#'     \item \code{"ln"} - Natural logarithm (default, most common in ecology)
#'     \item \code{"2"} - Log base 2 (information theory)
#'     \item \code{"10"} - Log base 10
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{Shannon} - Shannon diversity index values
#'   }
#'   Returns NA for stations with no individuals or all zero abundances.
#'
#' @details
#' The Shannon diversity index (H') is calculated as:
#'
#' \deqn{H' = -\sum_{i=1}^{S} p_i \log(p_i)}
#'
#' where:
#' \itemize{
#'   \item S = number of species
#'   \item p_i = proportion of individuals belonging to species i
#' }
#'
#' The index increases with both species richness and evenness. Higher values
#' indicate greater diversity. Typical values range from 0 (one species) to
#' 4-5 (very diverse communities).
#'
#' Zero abundances and NA values are automatically removed before calculation.
#'
#' @examples
#' # Example 1: Wide format (stations in rows, species in columns)
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 12, 3),
#'   Owenia = c(23, 18, 45),
#'   Nephtys = c(12, 34, 23),
#'   Lumbrineris = c(8, 15, 12)
#' )
#' Shannon(data_wide, base = "ln", shape = "w")
#'
#' # Example 2: Long format (species in rows, stations in columns)
#' data_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys", "Lumbrineris"),
#'   St1 = c(45, 23, 12, 8),
#'   St2 = c(12, 18, 34, 15),
#'   St3 = c(3, 45, 23, 12)
#' )
#' Shannon(data_long, base = "ln", shape = "l")
#'
#' # Example 3: Using different logarithm bases
#' Shannon(data_wide, base = "2", shape = "w")   # Information theory
#' Shannon(data_wide, base = "10", shape = "w")  # Common logarithm
#'
#' @references
#' Shannon, C.E. (1948). A mathematical theory of communication.
#' Bell System Technical Journal, 27(3), 379-423.
#'
#' Shannon, C.E. & Weaver, W. (1949). The Mathematical Theory of Communication.
#' University of Illinois Press.
#'
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' Gray, J.S. & Elliott, M. (2009). Ecology of Marine Sediments: From Science
#' to Management (2nd ed.). Oxford University Press, New York. 225 pp.
#'
#' @seealso
#' \code{\link{Simpson}} for Simpson's dominance index,
#' \code{\link{Pielou}} for Pielou's evenness,
#' \code{\link{Div}} for calculating multiple diversity indices at once
#'
#' @export
Shannon <- function(DF, base = "ln", shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!base %in% c("ln", "2", "10")) {
    stop("base must be 'ln', '2', or '10'")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Species <- DF[, 1]
    DF2 <- as.data.frame(t(DF[, -1]))
    colnames(DF2) <- Species
    rownames(DF2) <- colnames(DF[, -1])
    Stations <- rownames(DF2)

  } else {
    # Wide format: Stations in first column, species in other columns
    DF2 <- DF[, -1]
    Species <- colnames(DF2)
    Stations <- DF[, 1]
    rownames(DF2) <- Stations
  }

  # Calculate Shannon index for each station
  SV <- sapply(1:nrow(DF2), function(i) {
    x <- as.numeric(DF2[i, ])
    x <- x[!is.na(x) & x > 0]

    if (length(x) == 0) return(NA)

    pi <- x / sum(x)
    pi <- pi[pi > 0]

    if (base == "2") {
      -sum(pi * log2(pi))
    } else if (base == "10") {
      -sum(pi * log10(pi))
    } else {
      -sum(pi * log(pi))
    }
  })

  # Return results as data.frame
  Shannon_df <- data.frame(
    Stations = Stations,
    Shannon = SV,
    row.names = NULL
  )

  return(Shannon_df)
}


#' Calculate Species Richness
#'
#' Calculates the number of species (species richness) present in each station
#' or sampling unit for benthic community data.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{Richness} - Number of species present (abundance > 0)
#'   }
#'
#' @details
#' Species richness (S) is the simplest measure of biodiversity, representing
#' the total number of different species present in a sample. This function
#' converts abundance data to presence/absence and counts the number of species
#' with at least one individual.
#'
#' The function automatically handles:
#' \itemize{
#'   \item Zero abundances (not counted as present)
#'   \item NA values (ignored in the count)
#'   \item Negative values (returns an error)
#' }
#'
#' @examples
#' # Example 1: Wide format (stations in rows, species in columns)
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 0, 3),
#'   Owenia = c(23, 18, 0),
#'   Nephtys = c(0, 34, 23),
#'   Lumbrineris = c(8, 15, 12)
#' )
#' Richness(data_wide, shape = "w")
#'
#' # Example 2: Long format (species in rows, stations in columns)
#' data_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys", "Lumbrineris"),
#'   St1 = c(45, 23, 0, 8),
#'   St2 = c(0, 18, 34, 15),
#'   St3 = c(3, 0, 23, 12)
#' )
#' Richness(data_long, shape = "l")
#'
#' @references
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' Gray, J.S. & Elliott, M. (2009). Ecology of Marine Sediments: From Science
#' to Management (2nd ed.). Oxford University Press, New York. 225 pp.
#'
#' @export
Richness <- function(DF, shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Species <- DF[, 1]
    DF1 <- as.data.frame(t(DF[, -1]))
    colnames(DF1) <- Species
    rownames(DF1) <- colnames(DF[, -1])
    Stations <- rownames(DF1)

  } else {
    # Wide format: Stations in first column, species in other columns
    DF1 <- DF[, -1]
    Species <- colnames(DF1)
    Stations <- as.character(DF[, 1])
    rownames(DF1) <- Stations
  }

  # Validate non-negative abundances
  if (any(DF1 < 0, na.rm = TRUE)) {
    stop("Abundance values cannot be negative. Please check your data.")
  }

  # Convert to presence/absence
  PA <- DF1 > 0

  # Calculate richness (number of species present)
  richness <- rowSums(PA, na.rm = TRUE)

  # Return results as data.frame
  Richness_df <- data.frame(
    Stations = Stations,
    Richness = richness,
    row.names = NULL
  )

  return(Richness_df)
}


#' Calculate Maximum Shannon Diversity (Hmax)
#'
#' Calculates the maximum possible Shannon diversity (Hmax) for each station,
#' which represents the Shannon index value if all species were equally abundant.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{Hmax} - Maximum possible Shannon diversity (natural logarithm of species richness)
#'   }
#'
#' @details
#' Hmax represents the theoretical maximum value of Shannon diversity (H') that
#' would occur if all species in a sample were equally abundant. It is calculated as:
#'
#' \deqn{H_{max} = \ln(S)}
#'
#' where S is the number of species (species richness).
#'
#' This value is commonly used to calculate Pielou's evenness index (J'):
#'
#' \deqn{J' = H' / H_{max}}
#'
#' where H' is the observed Shannon diversity.
#'
#' @examples
#' # Example: Wide format
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 12, 3),
#'   Owenia = c(23, 18, 45),
#'   Nephtys = c(12, 34, 23),
#'   Lumbrineris = c(8, 15, 12)
#' )
#' Hmax(data_wide, shape = "w")
#'
#' # Compare with actual Shannon diversity
#' H <- Shannon(data_wide, shape = "w")
#' Hm <- Hmax(data_wide, shape = "w")
#' evenness <- H$Shannon / Hm$Hmax
#' data.frame(Station = H$Stations, H = H$Shannon, Hmax = Hm$Hmax, Evenness = evenness)
#'
#' @references
#' Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
#' consequences. Ecology, 54(2), 427-432.
#'
#' Pielou, E.C. (1966). The measurement of diversity in different types of
#' biological collections. Journal of Theoretical Biology, 13, 131-144.
#'
#' Shannon, C.E. & Weaver, W. (1949). The Mathematical Theory of Communication.
#' University of Illinois Press.
#'
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' @seealso \code{\link{Shannon}} for Shannon diversity index,
#'   \code{\link{Richness}} for species richness
#'
#' @export
Hmax <- function(DF, shape = "w") {

  # Calculate species richness using the Richness function
  S <- Richness(DF, shape = shape)

  # Calculate natural logarithm of richness
  LnS <- log(S$Richness)

  # Return results as data.frame
  Hmax_df <- data.frame(
    Stations = as.character(S$Stations),
    Hmax = LnS
  )

  return(Hmax_df)
}


#' Calculate Pielou's Evenness Index
#'
#' Calculates Pielou's evenness index (J'), which measures how evenly individuals
#' are distributed among species in a community.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param base Logarithm base for Shannon calculation. Options are:
#'   \itemize{
#'     \item \code{"ln"} - Natural logarithm (default, most common in ecology)
#'     \item \code{"2"} - Log base 2 (information theory)
#'     \item \code{"10"} - Log base 10
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{J} - Pielou's evenness index (values range from 0 to 1)
#'   }
#'
#' @details
#' Pielou's evenness index (J') measures how evenly individuals are distributed
#' among species. It is calculated as:
#'
#' \deqn{J' = H' / H_{max}}
#'
#' where:
#' \itemize{
#'   \item H' is the observed Shannon diversity index
#'   \item H_max is the maximum possible Shannon diversity (ln of species richness)
#' }
#'
#' Values range from 0 to 1:
#' \itemize{
#'   \item J' = 0: Maximum unevenness (one species dominates)
#'   \item J' = 1: Perfect evenness (all species equally abundant)
#' }
#'
#' Higher values indicate more even distribution of individuals among species.
#'
#' @examples
#' # Example 1: Community with high evenness
#' data_even <- data.frame(
#'   Station = c("St1"),
#'   Sp1 = 25, Sp2 = 25, Sp3 = 25, Sp4 = 25
#' )
#' Pielou(data_even, shape = "w")
#'
#' # Example 2: Community with low evenness (one dominant species)
#' data_uneven <- data.frame(
#'   Station = c("St2"),
#'   Sp1 = 97, Sp2 = 1, Sp3 = 1, Sp4 = 1
#' )
#' Pielou(data_uneven, shape = "w")
#'
#' # Example 3: Multiple stations
#' data_multi <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 90, 25),
#'   Owenia = c(23, 5, 25),
#'   Nephtys = c(12, 3, 25),
#'   Lumbrineris = c(8, 2, 25)
#' )
#' Pielou(data_multi, shape = "w")
#'
#' @references
#' Pielou, E.C. (1966). The measurement of diversity in different types of
#' biological collections. Journal of Theoretical Biology, 13, 131-144.
#' https://doi.org/10.1016/0022-5193(66)90013-0
#'
#' Pielou, E.C. (1975). Ecological Diversity. Wiley-Interscience, New York.
#'
#' Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
#' consequences. Ecology, 54(2), 427-432.
#'
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' @seealso \code{\link{Shannon}} for Shannon diversity index,
#'   \code{\link{Hmax}} for maximum Shannon diversity
#'
#' @export
Pielou <- function(DF, base = "ln", shape = "w") {

  # Calculate Shannon diversity
  H <- Shannon(DF, base = base, shape = shape)

  # Calculate maximum Shannon diversity
  HMax <- Hmax(DF, shape = shape)

  # Calculate Pielou's evenness
  JV <- H$Shannon / HMax$Hmax

  # Return results as data.frame
  J_df <- data.frame(
    Stations = H$Stations,
    J = JV
  )

  return(J_df)
}


#' Calculate Simpson's Dominance Index
#'
#' Calculates Simpson's dominance index (lambda), which measures the probability that
#' two randomly selected individuals belong to the same species.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{Dominance} - Simpson's dominance index (values range from 0 to 1)
#'   }
#'
#' @details
#' Simpson's dominance index (lambda) represents the probability that two individuals
#' randomly selected from a sample will belong to the same species. It is calculated as:
#'
#' \deqn{\lambda = \sum_{i=1}^{S} p_i^2}
#'
#' where:
#' \itemize{
#'   \item S = number of species
#'   \item p_i = proportion of individuals belonging to species i
#' }
#'
#' Values range from 0 to 1:
#' \itemize{
#'   \item lambda close to 0: Low dominance, high diversity (many species equally abundant)
#'   \item lambda close to 1: High dominance, low diversity (one species dominates)
#' }
#'
#' The Simpson diversity index (1 - lambda) can be calculated by subtracting this value from 1.
#'
#' @examples
#' # Example 1: High dominance (one species dominates)
#' data_dominated <- data.frame(
#'   Station = c("St1"),
#'   Sp1 = 90, Sp2 = 5, Sp3 = 3, Sp4 = 2
#' )
#' Dominance(data_dominated, shape = "w")  # High value (close to 1)
#'
#' # Example 2: Low dominance (species equally abundant)
#' data_even <- data.frame(
#'   Station = c("St2"),
#'   Sp1 = 25, Sp2 = 25, Sp3 = 25, Sp4 = 25
#' )
#' Dominance(data_even, shape = "w")  # Low value (0.25)
#'
#' # Example 3: Multiple stations
#' data_multi <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(80, 25, 10),
#'   Owenia = c(10, 25, 30),
#'   Nephtys = c(5, 25, 30),
#'   Lumbrineris = c(5, 25, 30)
#' )
#' Dominance(data_multi, shape = "w")
#'
#' @references
#' Simpson, E.H. (1949). Measurement of diversity. Nature, 163, 688.
#' https://doi.org/10.1038/163688a0
#'
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
#' consequences. Ecology, 54(2), 427-432.
#'
#' @seealso \code{\link{Shannon}} for Shannon diversity index,
#'   \code{\link{Pielou}} for evenness
#'
#' @export
Dominance <- function(DF, shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Species <- DF[, 1]
    DF3 <- as.data.frame(t(DF[, -1]))
    colnames(DF3) <- Species
    rownames(DF3) <- colnames(DF[, -1])
    Stations <- rownames(DF3)

  } else {
    # Wide format: Stations in first column, species in other columns
    DF3 <- DF[, -1]
    Species <- colnames(DF3)
    Stations <- as.character(DF[, 1])
    rownames(DF3) <- Stations
  }

  # Calculate proportions
  pi <- DF3 / rowSums(DF3, na.rm = TRUE)

  # Calculate Simpson's dominance index
  dominance <- rowSums(pi^2, na.rm = TRUE)

  # Return results as data.frame
  Dominance_df <- data.frame(
    Stations = Stations,
    Dominance = dominance
  )

  return(Dominance_df)
}


#' Calculate Simpson's Diversity Index
#'
#' Calculates Simpson's diversity index (1 - lambda), which measures the probability
#' that two randomly selected individuals belong to different species.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{Simpson} - Simpson's diversity index (values range from 0 to 1)
#'   }
#'
#' @details
#' Simpson's diversity index measures the probability that two individuals randomly
#' selected from a sample belong to different species. It is calculated as:
#'
#' \deqn{D = 1 - \lambda = 1 - \sum_{i=1}^{S} p_i^2}
#'
#' where:
#' \itemize{
#'   \item S = number of species
#'   \item p_i = proportion of individuals belonging to species i
#'   \item lambda = Simpson's dominance index
#' }
#'
#' Values range from 0 to 1:
#' \itemize{
#'   \item D close to 0: Low diversity (one species dominates)
#'   \item D close to 1: High diversity (many species equally abundant)
#' }
#'
#' This index is less sensitive to species richness and more sensitive to evenness
#' compared to Shannon's index. It gives more weight to abundant species.
#'
#' @examples
#' # Example 1: High diversity (species equally abundant)
#' data_diverse <- data.frame(
#'   Station = c("St1"),
#'   Sp1 = 25, Sp2 = 25, Sp3 = 25, Sp4 = 25
#' )
#' Simpson(data_diverse, shape = "w")  # High value (0.75)
#'
#' # Example 2: Low diversity (one species dominates)
#' data_dominated <- data.frame(
#'   Station = c("St2"),
#'   Sp1 = 90, Sp2 = 5, Sp3 = 3, Sp4 = 2
#' )
#' Simpson(data_dominated, shape = "w")  # Low value
#'
#' # Example 3: Multiple stations
#' data_multi <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(80, 25, 10),
#'   Owenia = c(10, 25, 30),
#'   Nephtys = c(5, 25, 30),
#'   Lumbrineris = c(5, 25, 30)
#' )
#' Simpson(data_multi, shape = "w")
#'
#' @references
#' Simpson, E.H. (1949). Measurement of diversity. Nature, 163, 688.
#' https://doi.org/10.1038/163688a0
#'
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
#' consequences. Ecology, 54(2), 427-432.
#'
#' @seealso \code{\link{Dominance}} for Simpson's dominance index (lambda),
#'   \code{\link{Shannon}} for Shannon diversity index
#'
#' @export
Simpson <- function(DF, shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Species <- DF[, 1]
    DF3 <- as.data.frame(t(DF[, -1]))
    colnames(DF3) <- Species
    rownames(DF3) <- colnames(DF[, -1])
    Stations <- rownames(DF3)

  } else {
    # Wide format: Stations in first column, species in other columns
    DF3 <- DF[, -1]
    Species <- colnames(DF3)
    Stations <- as.character(DF[, 1])
    rownames(DF3) <- Stations
  }

  # Calculate proportions
  pi <- DF3 / rowSums(DF3, na.rm = TRUE)

  # Calculate Simpson's diversity index (1 - dominance)
  simpson <- 1 - rowSums(pi^2, na.rm = TRUE)

  # Return results as data.frame
  Simpson_df <- data.frame(
    Stations = Stations,
    Simpson = simpson
  )

  return(Simpson_df)
}


#' Calculate Expected Species Richness (Rarefaction)
#'
#' Calculates the expected number of species in a standardized sample size
#' using rarefaction, allowing comparison of species richness across samples
#' with different sampling efforts.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param nsize Integer. Sample size for rarefaction (number of individuals to standardize to).
#'   Must be less than or equal to the total abundance in each station.
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{Esn} - Expected number of species for the standardized sample size (rounded to nearest integer)
#'   }
#'   Returns NA for stations where sample exceeds total abundance.
#'
#' @details
#' Rarefaction estimates the expected number of species in a random subsample
#' of individuals, standardizing richness estimates across samples with different
#' sample sizes or abundances. This method allows fair comparison of diversity
#' between samples with unequal sampling effort.
#'
#' The rarefaction method (Sanders 1968; Hurlbert 1971) calculates the expected
#' number of species (ESn) that would be observed if a smaller number (n) of
#' individuals were sampled from a community with known total species (S) and
#' total individuals (N). This provides a standardized richness measure (e.g., ES100
#' representing expected species in 100 individuals) that can be compared across
#' samples of varying sizes.
#'
#' For each species i, the probability of occurrence in a subsample of size n is
#' calculated as:
#'
#' \deqn{P_i = 1 - \frac{\binom{N - n_i}{n}}{\binom{N}{n}}}
#'
#' where:
#' \itemize{
#'   \item N = total number of individuals in the sample
#'   \item n = subsample size (sample parameter)
#'   \item n_i = number of individuals of species i
#' }
#'
#' The expected number of species (E(Sn)) is the sum of these probabilities
#' across all species. The calculation uses logarithmic transformations
#' (lchoose) to avoid numerical overflow with large samples.
#'
#' \strong{Important assumptions and limitations:}
#'
#' Rarefaction assumes that individuals are sampled independently of each other,
#' which corresponds to species being spatially randomly distributed (Poisson distribution).
#' In practice, this assumption is rarely met as most benthic species exhibit spatial
#' clustering or aggregation, which can be extreme. When spatial clustering occurs,
#' rarefaction tends to overestimate the expected number of species for smaller sample
#' sizes, making it a conservative measure of diversity.
#'
#' Despite this limitation, rarefaction remains a valuable tool for comparing diversity
#' across samples with different sampling efforts, particularly when combined with other
#' diversity metrics and ecological context.
#'
#' @examples
#' # Example 1: Rarefy to smallest nsize size
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 12, 3),
#'   Owenia = c(23, 8, 2),
#'   Nephtys = c(12, 5, 1),
#'   Lumbrineris = c(8, 3, 1)
#' )
#'
#' # Total abundances: St1=88, St2=28, St3=7
#' # Rarefy all to n=7 (smallest nsize)
#' Esn(data_wide, nsize = 7, shape = "w")
#'
#' # Example 2: Compare richness at different rarefaction levels
#' Esn(data_wide, nsize = 10, shape = "w")
#' Esn(data_wide, nsize = 20, shape = "w")
#'
#' # Example 3: Long format
#' data_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys"),
#'   St1 = c(45, 23, 12),
#'   St2 = c(12, 8, 5)
#' )
#' Esn(data_long, nsize = 15, shape = "l")
#'
#' @references
#' Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique and
#' alternative parameters. Ecology, 52(4), 577-586.
#' https://doi.org/10.2307/1934145
#'
#' Sanders, H.L. (1968). Marine benthic diversity: a comparative study.
#' The American Naturalist, 102(925), 243-282.
#'
#' Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity: procedures
#' and pitfalls in the measurement and comparison of species richness.
#' Ecology Letters, 4(4), 379-391.
#'
#' Heck, K.L., van Belle, G. & Simberloff, D. (1975). Explicit calculation of
#' the rarefaction diversity measurement and the determination of sufficient
#' sample size. Ecology, 56(6), 1459-1461.
#'
#' Clarke, K.R., Gorley, R.N., Somerfield, P.J. & Warwick, R.M. (2014).
#' Change in Marine Communities: An Approach to Statistical Analysis and
#' Interpretation (3rd ed.). PRIMER-E, Plymouth.
#'
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' @export
Esn <- function(DF, nsize, shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!is.numeric(nsize) || length(nsize) != 1 || nsize <= 0) {
    stop("nsize must be a single positive number")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    Species <- DF[, 1]
    DF4 <- as.data.frame(t(DF[, -1]))
    colnames(DF4) <- Species
    rownames(DF4) <- colnames(DF[, -1])
    Stations <- rownames(DF4)
  } else {
    DF4 <- DF[, -1]
    Species <- colnames(DF4)
    Stations <- as.character(DF[, 1])
    rownames(DF4) <- Stations
  }

  # Calculate total abundance per station
  N <- rowSums(DF4, na.rm = TRUE)

  # Check if nsize exceeds any station's abundance
  if (any(nsize > N, na.rm = TRUE)) {
    warning(paste("nsize =", nsize, "exceeds total abundance for some stations.\nThese will return NA."))
  }

  # Initialize result vector
  Esn_val <- numeric(nrow(DF4))

  # Calculate rarefaction for each station
  for (j in 1:nrow(DF4)) {

    # Skip if nsize exceeds station abundance or abundance is zero
    if (nsize > N[j] || N[j] == 0) {
      Esn_val[j] <- NA
      next
    }

    # Initialize presence probability vector
    presence <- numeric(ncol(DF4))

    # Calculate probability for each species
    for (i in 1:ncol(DF4)) {

      ni <- DF4[j, i]

      # Skip absent or NA species
      if (is.na(ni) || ni == 0) {
        presence[i] <- 0
        next
      }

      # If species is very abundant, it will definitely be present
      if (N[j] - ni < nsize) {
        presence[i] <- 1
      } else {
        # Calculate probability of absence using combinatorics
        # P(absence) = C(N-ni, n) / C(N, n)
        # Using log-scale to avoid numerical overflow
        absence <- exp(lchoose(N[j] - ni, nsize) - lchoose(N[j], nsize))

        # Handle numerical issues
        if (is.nan(absence) || is.infinite(absence)) {
          absence <- 0
        }

        # P(presence) = 1 - P(absence)
        presence[i] <- 1 - absence
      }
    }

    # Sum probabilities to get expected richness
    Esn_val[j] <- sum(presence)
  }

  # Return results as data.frame
  Esn_df <- data.frame(
    Stations = Stations,
    Esn = round(Esn_val),
    row.names = NULL
  )

  return(Esn_df)
}


#' Beta Diversity Partitioning (Baselga Framework)
#'
#' Calculates multiple-site beta diversity partitioning into turnover and
#' nestedness (or balanced variation and gradient) components using the
#' framework of Baselga (2010, 2012, 2013). Returns a summary table with
#' overall partitioning values across all sites.
#'
#' @param DF A data.frame containing species data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param index Type of beta diversity index to use:
#'   \itemize{
#'     \item \code{"sorensen"} - Sorensen dissimilarity (presence-absence, default)
#'     \item \code{"jaccard"} - Jaccard dissimilarity (presence-absence)
#'     \item \code{"bray"} - Bray-Curtis dissimilarity (abundance-based)
#'     \item \code{"ruzicka"} - Ruzicka dissimilarity (abundance-based)
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return A data.frame with three rows and three columns:
#'   \itemize{
#'     \item \code{Component} - Name of the beta diversity component
#'     \item \code{Beta} - Numeric value of the component (0 to 1)
#'     \item \code{Interpretation} - Brief explanation of what the component measures
#'   }
#'
#'   The three components depend on the index used:
#'
#'   For Sorensen: Total (betaSOR), Turnover (betaSIM), Nestedness (betaSNE)
#'
#'   For Jaccard: Total (betaJAC), Turnover (betaJTU), Nestedness (betaJNE)
#'
#'   For Bray-Curtis: Total (betaBRA), Balanced (betaBAL), Gradient (betaGRA)
#'
#'   For Ruzicka: Total (betaRUZ), Balanced (betaBAL), Gradient (betaGRA)
#'
#' @details
#' This function calculates multiple-site beta diversity partitioning following
#' the additive framework of Baselga (2010, 2012, 2013). It returns a single set
#' of values representing the overall partitioning across all sites in the dataset.
#'
#' \strong{For presence-absence data (Sorensen/Jaccard):}
#'
#' Total beta diversity is additively partitioned as:
#' \deqn{\beta_{total} = \beta_{turnover} + \beta_{nestedness}}
#'
#' \itemize{
#'   \item \strong{Total (betaSOR/betaJAC)} - Overall compositional dissimilarity among all sites.
#'     Values close to 1 indicate high differentiation; close to 0 indicate similar communities.
#'   \item \strong{Turnover (betaSIM/betaJTU)} - Measures species replacement between sites.
#'     High values indicate communities differ primarily due to species substitution.
#'     Reflects deterministic processes like environmental filtering, dispersal limitation,
#'     or competitive exclusion.
#'   \item \strong{Nestedness (betaSNE/betaJNE)} - Measures species loss/gain without replacement.
#'     High values indicate some sites contain subsets of species from other sites.
#'     Reflects stochastic processes like colonization-extinction dynamics, habitat loss,
#'     or ordered extinctions along environmental gradients.
#' }
#'
#' \strong{For abundance data (Bray-Curtis/Ruzicka):}
#'
#' Total beta diversity is additively partitioned as:
#' \deqn{\beta_{total} = \beta_{balanced} + \beta_{gradient}}
#'
#' \itemize{
#'   \item \strong{Total (betaBRA/betaRUZ)} - Overall abundance-based dissimilarity among all sites.
#'   \item \strong{Balanced variation (betaBAL)} - Abundance changes balanced across species.
#'     Analogous to turnover in presence-absence framework. Reflects reciprocal abundance
#'     changes among species (some increase while others decrease).
#'   \item \strong{Abundance gradient (betaGRA)} - Unidirectional abundance changes.
#'     Analogous to nestedness in presence-absence framework. Reflects consistent
#'     abundance gradients across sites (e.g., productivity or disturbance gradients).
#' }
#'
#' \strong{Interpretation guidelines:}
#' \itemize{
#'   \item All values range from 0 to 1
#'   \item betatotal = betaturnover + betanestedness (additive partitioning)
#'   \item If Turnover/Balanced >> Nestedness/Gradient: Replacement processes dominate
#'   \item If Nestedness/Gradient >> Turnover/Balanced: Nested patterns or gradients dominate
#'   \item Balanced partitioning suggests both processes contribute equally
#' }
#'
#' \strong{Choosing an index:}
#' \itemize{
#'   \item Use \strong{Sorensen} for presence-absence data (most common, insensitive to richness differences)
#'   \item Use \strong{Jaccard} for presence-absence when richness differences matter
#'   \item Use \strong{Bray-Curtis} for abundance data (most common in ecology)
#'   \item Use \strong{Ruzicka} for abundance data (quantitative form of Jaccard)
#' }
#'
#' Understanding the partitioning helps identify ecological processes structuring
#' benthic communities and guides interpretation of diversity patterns.
#'
#' @examples
#' # Example 1: Presence-absence with Sorensen (default)
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3", "St4"),
#'   Capitella = c(1, 1, 0, 0),
#'   Owenia = c(1, 0, 1, 0),
#'   Nephtys = c(0, 0, 1, 1),
#'   Lumbrineris = c(0, 1, 1, 0)
#' )
#' BetaPart(data_wide, index = "sorensen", shape = "w")
#'
#' # Example 2: Presence-absence with Jaccard
#' result_jac <- BetaPart(data_wide, index = "jaccard", shape = "w")
#' print(result_jac)
#'
#' # Example 3: Abundance-based with Bray-Curtis
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3", "St4"),
#'   Capitella = c(45, 12, 0, 0),
#'   Owenia = c(23, 18, 45, 50),
#'   Nephtys = c(12, 0, 23, 10),
#'   Lumbrineris = c(0, 15, 12, 5)
#' )
#' result_bray <- BetaPart(data_wide, index = "bray", shape = "w")
#' print(result_bray)
#'
#' # Example 4: Abundance-based with Ruzicka
#' BetaPart(data_wide, index = "ruzicka", shape = "w")
#'
#' # Example 5: Long format data
#' data_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys"),
#'   St1 = c(1, 0, 1),
#'   St2 = c(1, 1, 0),
#'   St3 = c(0, 0, 1)
#' )
#' BetaPart(data_long, index = "sorensen", shape = "l")
#'
#' # Example 6: Interpreting results
#' result <- BetaPart(data_wide, index = "bray", shape = "w")
#' # High Balanced (betaBAL) suggests reciprocal abundance changes
#' # High Gradient (betaGRA) suggests unidirectional abundance patterns
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components of
#' beta diversity. Global Ecology and Biogeography, 19(1), 134-143.
#' https://doi.org/10.1111/j.1466-8238.2009.00490.x
#'
#' Baselga, A. (2012). The relationship between species replacement, dissimilarity
#' derived from nestedness, and nestedness. Global Ecology and Biogeography,
#' 21(12), 1223-1232. https://doi.org/10.1111/j.1466-8238.2011.00756.x
#'
#' Baselga, A. (2013). Separating the two components of abundance-based
#' dissimilarity: balanced changes in abundance vs. abundance gradients.
#' Methods in Ecology and Evolution, 4(6), 552-557.
#' https://doi.org/10.1111/2041-210X.12029
#'
#' Baselga, A. & Orme, C.D.L. (2012). betapart: an R package for the study
#' of beta diversity. Methods in Ecology and Evolution, 3(5), 808-812.
#' https://doi.org/10.1111/2041-210X.12110
#'
#' @seealso \code{\link{BetaDiv}} for LCBD/SCBD analysis using adespatial,
#'   \code{\link[betapart]{beta.multi}} for presence-absence partitioning,
#'   \code{\link[betapart]{beta.multi.abund}} for abundance-based partitioning
#'
#' @export
BetaPart <- function(DF, index = "sorensen", shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!index %in% c("sorensen", "jaccard", "bray", "ruzicka")) {
    stop("index must be 'sorensen', 'jaccard', 'bray', or 'ruzicka'")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Species <- DF[, 1]
    DF6 <- as.data.frame(t(DF[, -1]))
    colnames(DF6) <- Species
    rownames(DF6) <- colnames(DF[, -1])
  } else {
    # Wide format: Stations in first column, species in other columns
    DF6 <- DF[, -1]
    rownames(DF6) <- DF[, 1]
  }

  # Calculate beta diversity partitioning based on index
  if (index == "sorensen") {
    # Presence-absence: Sorensen
    result <- betapart::beta.multi(DF6, index.family = "sorensen")

    # Create summary
    summary_df <- data.frame(
      Component = c("Total (BetaSOR)", "Turnover (BetaSIM)", "Nestedness (BetaSNE)"),
      Beta = c(
        result$beta.SOR,
        result$beta.SIM,
        result$beta.SNE
      ),
      Interpretation = c(
        "Overall dissimilarity between sites",
        "Dissimilarity due to species replacement",
        "Dissimilarity due to species loss/gain"
      )
    )

  } else if (index == "jaccard") {
    # Presence-absence: Jaccard
    result <- betapart::beta.multi(DF6, index.family = "jaccard")

    # Create summary
    summary_df <- data.frame(
      Component = c("Total (BetaJAC)", "Turnover (BetaJTU)", "Nestedness (BetaJNE)"),
      Beta = c(
        result$beta.JAC,
        result$beta.JTU,
        result$beta.JNE
      ),
      Interpretation = c(
        "Overall dissimilarity between sites",
        "Dissimilarity due to species replacement",
        "Dissimilarity due to species loss/gain"
      )
    )

  } else if (index == "bray") {
    # Abundance: Bray-Curtis
    result <- betapart::beta.multi.abund(DF6, index.family = "bray")

    # Create summary
    summary_df <- data.frame(
      Component = c("Total (BetaBRA)", "Balanced (BetaBAL)", "Gradient (BetaGRA)"),
      Beta = c(
        result$beta.BRAY,
        result$beta.BRAY.BAL,
        result$beta.BRAY.GRA
      ),
      Interpretation = c(
        "Overall dissimilarity between sites",
        "Dissimilarity due to balanced abundance turnover",
        "Dissimilarity due to abundance gradients"
      )
    )

  } else if (index == "ruzicka") {
    # Abundance: Ruzicka
    result <- betapart::beta.multi.abund(DF6, index.family = "ruzicka")

    # Create summary
    summary_df <- data.frame(
      Component = c("Total (BetaRUZ)", "Balanced (BetaBAL)", "Gradient (BetaGRA)"),
      Beta = c(
        result$beta.RUZ,
        result$beta.RUZ.BAL,
        result$beta.RUZ.GRA
      ),
      Interpretation = c(
        "Overall dissimilarity between sites",
        "Dissimilarity due to balanced abundance turnover",
        "Dissimilarity due to abundance gradients"
      )
    )
  }

  return(summary_df)
}


#' Calculate Beta Diversity Indices (LCBD and SCBD)
#'
#' Calculates Local Contribution to Beta Diversity (LCBD) for stations and/or
#' Species Contribution to Beta Diversity (SCBD) for species, which measure
#' the uniqueness of sites and the contribution of species to beta diversity.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param output Character. Type of output to return:
#'   \itemize{
#'     \item \code{"LCBD"} - Returns only LCBD indices for stations (default)
#'     \item \code{"SCBD"} - Returns only SCBD indices for species
#'     \item \code{"both"} - Returns a list containing both LCBD and SCBD
#'   }
#' @param method Dissimilarity method for beta diversity partitioning. Default is "hellinger".
#'   Available transformations include:
#'   \itemize{
#'     \item \code{"hellinger"} - Hellinger transformation (default). Recommended for abundance data.
#'       Gives appropriate weight to rare species. Formula: \eqn{\sqrt{p_{ij}}} where p is proportion.
#'     \item \code{"chord"} - Chord transformation. Suitable for compositional data.
#'       Normalizes by Euclidean norm.
#'     \item \code{"profile"} - Chi-square transformation. Emphasizes relative abundances.
#'       Suitable for presence-absence derived data.
#'     \item \code{"chisquare"} - Chi-square distance. For count data with varying totals.
#'     \item \code{"total"} - Total abundance transformation. No transformation applied.
#'     \item \code{"max"} - Maximum transformation. Divides by maximum value.
#'     \item \code{"pa"} - Presence-absence transformation. Converts to binary (0/1).
#'     \item \code{"Wisconsin"} - Wisconsin double standardization. Species maxima then site totals.
#'     \item \code{"log.chord"} - Log-chord transformation. Log then chord.
#'     \item \code{"percentdiff"} - Percentage difference (Bray-Curtis). For abundance data.
#'     \item \code{"ruzicka"} - Ruzicka index. Quantitative form of Jaccard.
#'   }
#'   For benthic abundance data, \strong{"hellinger"} is recommended as it appropriately
#'   balances the influence of rare and abundant species.
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#'
#' @return Depending on \code{output}:
#'   \itemize{
#'     \item If \code{output = "LCBD"}: A data.frame with columns:
#'       \itemize{
#'         \item \code{Stations} - Station or sample names
#'         \item \code{LCBD} - Local Contribution to Beta Diversity index
#'         \item \code{pLCBD} - P-value from permutation test (999 permutations)
#'         \item \code{adpLCBD} - Adjusted p-value (Holm correction)
#'       }
#'     \item If \code{output = "SCBD"}: A data.frame with columns:
#'       \itemize{
#'         \item \code{Species} - Species names
#'         \item \code{SCBD} - Species Contribution to Beta Diversity index
#'       }
#'     \item If \code{output = "both"}: A list with two elements:
#'       \itemize{
#'         \item \code{$LCBD} - Data.frame with LCBD values
#'         \item \code{$SCBD} - Data.frame with SCBD values
#'       }
#'   }
#'
#' @details
#' This function calculates beta diversity partitioning using the framework of
#' Legendre & De Caceres (2013), providing two complementary perspectives:
#'
#' \strong{LCBD (Local Contribution to Beta Diversity):}
#'
#' LCBD indices measure the degree of uniqueness of each sampling unit (station)
#' in terms of community composition. High LCBD values indicate sites that are
#' compositionally different from the average, which may represent:
#' \itemize{
#'   \item Ecotones or transition zones
#'   \item Disturbed or impacted sites
#'   \item Sites with unique environmental conditions
#'   \item Biodiversity hotspots
#' }
#'
#' The sum of all LCBD values equals the total beta diversity (BDtotal).
#' Statistical significance is assessed through permutation tests (999 permutations).
#' Adjusted p-values account for multiple testing using the Holm correction.
#'
#' \strong{SCBD (Species Contribution to Beta Diversity):}
#'
#' SCBD indices measure how much each species contributes to overall beta diversity.
#' High SCBD values indicate species that:
#' \itemize{
#'   \item Vary greatly in abundance across sites
#'   \item Are good indicators of environmental gradients
#'   \item Differentiate community types
#'   \item May be keystone or indicator species
#' }
#'
#' Species with high SCBD are responsible for differentiating sites, while species
#' with low SCBD are more uniformly distributed. The sum of all SCBD values equals 1.
#'
#' \strong{Interpretation together:}
#'
#' Using both indices provides a complete picture:
#' \itemize{
#'   \item LCBD identifies \emph{which sites} are different
#'   \item SCBD identifies \emph{which species} cause those differences
#' }
#'
#' This combined approach is particularly useful for understanding community patterns,
#' identifying indicator species, and detecting environmental impacts.
#'
#' This function is a wrapper around \code{adespatial::beta.div()}, adapted for
#' flexible data input formats commonly used in benthic ecology.
#'
#' @examples
#' # Example data
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3", "St4"),
#'   Capitella = c(45, 12, 3, 0),
#'   Owenia = c(23, 18, 45, 50),
#'   Nephtys = c(12, 34, 23, 10),
#'   Lumbrineris = c(8, 15, 12, 5)
#' )
#'
#' # Example 1: Get only LCBD (station uniqueness)
#' lcbd_results <- BetaDiv(data_wide, output = "LCBD")
#' lcbd_results
#'
#' # Example 2: Get only SCBD (species contribution)
#' scbd_results <- BetaDiv(data_wide, output = "SCBD")
#' scbd_results
#'
#' # Example 3: Get both LCBD and SCBD
#' results <- BetaDiv(data_wide, output = "both")
#' results$LCBD  # Station indices
#' results$SCBD  # Species indices
#'
#' # Example 4: Using different transformation
#' BetaDiv(data_wide, output = "both", method = "chord")
#'
#'
#' # Example 5: Long format data
#' data_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys"),
#'   St1 = c(45, 23, 12),
#'   St2 = c(12, 18, 34),
#'   St3 = c(3, 45, 23)
#' )
#' BetaDiv(data_long, output = "both", shape = "l")
#'
#' @references
#' Legendre, P. & De Caceres, M. (2013). Beta diversity as the variance of
#' community data: dissimilarity coefficients and partitioning. Ecology Letters,
#' 16(8), 951-963. https://doi.org/10.1111/ele.12141
#'
#' Legendre, P. (2014). Interpreting the replacement and richness difference
#' components of beta diversity. Global Ecology and Biogeography, 23(11), 1324-1334.
#'
#' Borcard, D., Gillet, F. & Legendre, P. (2018). Numerical Ecology with R (2nd ed.).
#' Springer.
#'
#' @seealso \code{\link[adespatial]{beta.div}} for the underlying function,
#'   \code{\link{Shannon}} for alpha diversity,
#'   \code{\link{Richness}} for species richness
#'
#' @export
BetaDiv <- function(DF, output = "LCBD", method = "hellinger", shape = "w"){

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!output %in% c("LCBD", "SCBD", "both")) {
    stop("output must be 'LCBD', 'SCBD', or 'both'")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Species <- DF[, 1]
    DF5 <- as.data.frame(t(DF[, -1]))
    colnames(DF5) <- Species
    rownames(DF5) <- colnames(DF[, -1])
    Stations <- rownames(DF5)
  } else {
    # Wide format: Stations in first column, species in other columns
    DF5 <- DF[, -1]
    Species <- colnames(DF5)
    Stations <- as.character(DF[, 1])
    rownames(DF5) <- Stations
  }

  # Calculate beta diversity partitioning
  results <- adespatial::beta.div(Y = DF5, method = method, nperm = 999)

  # Prepare LCBD data.frame
  LCBD_df <- data.frame(
    Stations = Stations,
    LCBD = results$LCBD,
    pLCBD = results$p.LCBD,
    adpLCBD = results$p.adj,
    row.names = NULL
  )

  # Prepare SCBD data.frame
  SCBD_df <- data.frame(
    Species = names(results$SCBD),
    SCBD = as.numeric(results$SCBD),
    row.names = NULL
  )

  # Return based on output parameter
  if (output == "LCBD") {
    return(LCBD_df)
  } else if (output == "SCBD") {
    return(SCBD_df)
  } else {
    # Return both as a list
    return(list(
      LCBD = LCBD_df,
      SCBD = SCBD_df
    ))
  }
}


#' Calculate Multiple Diversity Indices
#'
#' Calculates a comprehensive set of diversity indices for benthic community data
#' in a single function call, including alpha diversity, rarefaction, and beta
#' diversity components.
#'
#' @param DF A data.frame containing species abundance data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param base Logarithm base for Shannon calculation. Options are:
#'   \itemize{
#'     \item \code{"ln"} - Natural logarithm (default, most common in ecology)
#'     \item \code{"2"} - Log base 2 (information theory)
#'     \item \code{"10"} - Log base 10
#'   }
#'   See \code{\link{Shannon}} for details.
#' @param n Integer. Sample size for rarefaction. Must be less than or equal
#'   to the smallest total abundance across stations. See \code{\link{Esn}} for details.
#' @param method Dissimilarity method for beta diversity (LCBD) calculation.
#'   Default is "hellinger". See \code{\link{BetaDiv}} for available methods.
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns (first column = species names)
#'   }
#' @param simpson Logical. Should Simpson's diversity index (1-lambda) be included?
#'   Default is FALSE. Set to TRUE to include Simpson diversity in the output.
#'
#' @return A data.frame with one row per station and the following columns:
#'   \itemize{
#'     \item \code{Stations} - Station or sample names
#'     \item \code{N} - Total abundance (number of individuals)
#'     \item \code{Richness} - Number of species (S)
#'     \item \code{Esn} - Expected species richness at standardized sample size (rarefaction)
#'     \item \code{Shannon} - Shannon diversity index (H')
#'     \item \code{Hmax} - Maximum possible Shannon diversity
#'     \item \code{J} - Pielou's evenness index (J')
#'     \item \code{Dominance} - Simpson's dominance index (lambda)
#'     \item \code{Simpson} - Simpson's diversity index (1-lambda) (only if simpson = TRUE)
#'     \item \code{LCBD} - Local contribution to beta diversity
#'     \item \code{pLCBD} - P-value for LCBD (from permutation test)
#'     \item \code{adpLCBD} - Adjusted p-value for LCBD (Holm correction)
#'   }
#'
#' @details
#' This function is a convenient wrapper that calculates multiple diversity indices
#' in a single call, providing a comprehensive characterization of benthic community
#' diversity at each station. It combines:
#'
#' \strong{Alpha diversity metrics:}
#' \itemize{
#'   \item \strong{Abundance (N)} - Total number of individuals
#'   \item \strong{Richness (S)} - Number of species present
#'   \item \strong{Rarefied richness (Esn)} - Expected species richness standardized
#'     to sample size n, allowing fair comparison across samples with different
#'     sampling efforts
#'   \item \strong{Shannon diversity (H')} - Accounts for both richness and evenness,
#'     incorporating the proportional abundance of each species
#'   \item \strong{Maximum diversity (Hmax)} - Theoretical maximum Shannon diversity
#'     if all species were equally abundant
#'   \item \strong{Evenness (J')} - Pielou's index measuring how evenly individuals
#'     are distributed among species (H'/Hmax)
#'   \item \strong{Dominance (lambda)} - Simpson's dominance index, the probability that
#'     two randomly selected individuals belong to the same species
#'   \item \strong{Simpson (1-lambda)} - Simpson's diversity index, the probability that
#'     two randomly selected individuals belong to different species (optional)
#' }
#'
#' \strong{Beta diversity metric:}
#' \itemize{
#'   \item \strong{LCBD} - Local Contribution to Beta Diversity, measuring how
#'     compositionally unique each station is compared to the overall dataset
#' }
#'
#' All individual functions are called internally with the specified parameters.
#' Error handling ensures that if one index fails to calculate, others will still
#' be returned with NA values for the failed index.
#'
#' This comprehensive approach is particularly useful for:
#' \itemize{
#'   \item Routine diversity assessments
#'   \item Comparative studies across multiple stations
#'   \item Environmental impact assessments
#'   \item Monitoring programs
#'   \item Ecological status evaluations
#' }
#'
#' @examples
#' # Example 1: Calculate all diversity indices (default, no Simpson)
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 12, 3),
#'   Owenia = c(23, 18, 45),
#'   Nephtys = c(12, 34, 23),
#'   Lumbrineris = c(8, 15, 12)
#' )
#'
#' # Calculate all indices (rarefy to n=60, natural log for Shannon)
#' diversity_results <- Div(data_wide, base = "ln", n = 60,
#'                          method = "hellinger", shape = "w")
#' diversity_results
#'
#' # Example 2: Include Simpson diversity index
#' Div(data_wide, base = "ln", n = 60, method = "hellinger",
#'     shape = "w", simpson = TRUE)
#'
#' # Example 3: Using log base 2 for Shannon
#' Div(data_wide, base = "2", n = 60, method = "hellinger", shape = "w")
#'
#' # Example 4: Using chord transformation for beta diversity
#' Div(data_wide, base = "ln", n = 60, method = "chord", shape = "w")
#'
#' # Example 5: Long format data
#' data_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys"),
#'   St1 = c(45, 23, 12),
#'   St2 = c(12, 18, 34),
#'   St3 = c(3, 45, 23)
#' )
#' Div(data_long, base = "ln", n = 50, method = "hellinger", shape = "l")
#'
#' @references
#' Magurran, A.E. (2004). Measuring Biological Diversity. Blackwell Publishing.
#'
#' Gray, J.S. & Elliott, M. (2009). Ecology of Marine Sediments: From Science
#' to Management (2nd ed.). Oxford University Press, New York. 225 pp.
#'
#' Legendre, P. & De Caceres, M. (2013). Beta diversity as the variance of
#' community data: dissimilarity coefficients and partitioning. Ecology Letters,
#' 16(8), 951-963. https://doi.org/10.1111/ele.12141
#'
#' Shannon, C.E. & Weaver, W. (1949). The Mathematical Theory of Communication.
#' University of Illinois Press.
#'
#' Simpson, E.H. (1949). Measurement of diversity. Nature, 163, 688.
#' https://doi.org/10.1038/163688a0
#'
#' @seealso \code{\link{Shannon}}, \code{\link{Richness}}, \code{\link{Esn}},
#'   \code{\link{Pielou}}, \code{\link{Dominance}}, \code{\link{Simpson}},
#'   \code{\link{BetaDiv}}
#'
#' @export
Div <- function(DF, base = "ln", n, method = "hellinger", shape = "w", simpson = FALSE) {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (missing(n)) {
    stop("n (sample size for rarefaction) is required")
  }

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("n must be a single positive number")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  if (!is.logical(simpson) || length(simpson) != 1) {
    stop("simpson must be TRUE or FALSE")
  }

  # Data preparation based on format
  if (shape == "l") {
    Species <- DF[, 1]
    DF5 <- as.data.frame(t(DF[, -1]))
    colnames(DF5) <- Species
    rownames(DF5) <- colnames(DF[, -1])
    Stations <- rownames(DF5)
  } else {
    DF5 <- DF[, -1]
    Species <- colnames(DF5)
    Stations <- DF[, 1]
    rownames(DF5) <- Stations
  }

  # Calculate total abundance (always succeeds)
  N <- rowSums(DF5)
  N_df <- data.frame(Stations = Stations, N = N)

  # Calculate diversity indices with error handling

  # Richness
  S <- tryCatch(
    Richness(DF, shape = shape),
    error = function(e) {
      warning("Richness calculation failed: ", e$message)
      data.frame(Stations = Stations, Richness = NA)
    }
  )

  # Rarefaction
  EsnR <- tryCatch(
    Esn(DF, nsize = n, shape = shape),
    error = function(e) {
      warning("Rarefaction (Esn) calculation failed: ", e$message)
      data.frame(Stations = Stations, Esn = NA)
    }
  )

  # Shannon
  H <- tryCatch(
    Shannon(DF, base = base, shape = shape),
    error = function(e) {
      warning("Shannon calculation failed: ", e$message)
      data.frame(Stations = Stations, Shannon = NA)
    }
  )

  # Hmax
  HMax <- tryCatch(
    Hmax(DF, shape = shape),
    error = function(e) {
      warning("Hmax calculation failed: ", e$message)
      data.frame(Stations = Stations, Hmax = NA)
    }
  )

  # Pielou
  J <- tryCatch(
    Pielou(DF, base = base, shape = shape),
    error = function(e) {
      warning("Pielou calculation failed: ", e$message)
      data.frame(Stations = Stations, J = NA)
    }
  )

  # Dominance
  L <- tryCatch(
    Dominance(DF, shape = shape),
    error = function(e) {
      warning("Dominance calculation failed: ", e$message)
      data.frame(Stations = Stations, Dominance = NA)
    }
  )

  # LCBD
  LCBD <- tryCatch(
    BetaDiv(DF, output = "LCBD", method = method, shape = shape),
    error = function(e) {
      warning("LCBD calculation failed: ", e$message)
      data.frame(Stations = Stations, LCBD = NA, pLCBD = NA, adpLCBD = NA)
    }
  )

  # Build list of results (always include base indices)
  Div_list <- list(N_df, S, EsnR, H, HMax, J, L, LCBD)

  # Optionally add Simpson
  if (simpson) {
    Simp <- tryCatch(
      Simpson(DF, shape = shape),
      error = function(e) {
        warning("Simpson calculation failed: ", e$message)
        data.frame(Stations = Stations, Simpson = NA)
      }
    )
    # Insert Simpson after Dominance (before LCBD)
    Div_list <- list(N_df, S, EsnR, H, HMax, J, L, Simp, LCBD)
  }

  # Combine all results
  Div_df <- Reduce(function(x, y) merge(x, y, by = "Stations", all.x = TRUE, sort = FALSE),
                   Div_list)

  return(Div_df)
}

