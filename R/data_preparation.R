#' Transform Data to Achieve Normality
#'
#' Applies appropriate transformation methods based on data type (environmental,
#' biological, percentage/coverage, or presence-absence) to achieve normality
#' for parametric statistical analyses.
#'
#' @param DF A data.frame containing variables to transform. The first column
#'   should contain sample/station names, remaining columns are numeric variables.
#' @param data_type Character. Type of data to determine appropriate transformations:
#'   \itemize{
#'     \item \code{"env"} - Environmental variables (default): sqrt → log → fourth root → chord → scale
#'     \item \code{"bio"} - Biological/abundance data: sqrt → fourth root → log → chord → scale
#'     \item \code{"percent"} - Percentage/proportion data in the range of 0 to 1 or in the range of 0 to 100: arcsine-sqrt only
#'     \item \code{"pa"} - Presence-absence data: converts to binary (0/1), no normality testing
#'   }
#' @param alpha Significance level for normality test (default = 0.05).
#' @param min_sd Minimum standard deviation threshold. Variables with SD below
#'   this value are excluded (default = 0, keeps all variables with any variation).
#' @param test Normality test to use: "shapiro" (default) or "ks".
#'
#' @return A list with six elements:
#'   \itemize{
#'     \item \code{normal_data} - Data frame with all normalized variables (ready
#'       for parametric analysis)
#'     \item \code{non_normal_data} - Data frame with variables that could not be
#'       normalized (use non-parametric methods)
#'     \item \code{transformation_report} - Detailed report showing which transformation
#'       was applied to each variable
#'     \item \code{summary} - Summary statistics: total variables, normalized,
#'       non-normalized, and method counts
#'     \item \code{normal_vars} - Character vector of normalized variable names
#'     \item \code{non_normal_vars} - Character vector of non-normalized variable names
#'   }
#'
#' @details
#' This function applies data-type-specific transformation workflows to achieve
#' normality for parametric statistical analyses.
#'
#' \strong{Transformation Sequences by Data Type:}
#'
#' \describe{
#'   \item{\strong{Environmental ("env"):}}{
#'     Sequence: Initial test → Square root → Log(x+1) → Fourth root → Chord → Z-score
#'
#'     \strong{Rationale:} Environmental variables (temperature, salinity, nutrients,
#'     depth, etc.) typically show moderate skewness and benefit from progressive
#'     transformations. Start with mild (sqrt), progress to moderate (log), then
#'     aggressive (fourth root), followed by standardization methods.
#'
#'     \strong{Typical variables:} Temperature, salinity, depth, pH, dissolved oxygen,
#'     nutrients (NO3, PO4), chlorophyll, turbidity, grain size parameters
#'   }
#'
#'   \item{\strong{Biological/Abundance ("bio"):}}{
#'     Sequence: Initial test → Square root → Fourth root → Log(x+1) → Chord → Z-score
#'
#'     \strong{Rationale:} Biological abundance data are typically highly right-skewed
#'     with many zeros and few very abundant species. Fourth root is preferred in
#'     benthic ecology as it down-weights dominant species while retaining
#'     community structure better than log transformation.
#'
#'     \strong{Typical variables:} Species abundances, densities (ind/m²), biomass,
#'     counts per taxonomic group
#'
#'     \strong{Note:} For species matrices intended for beta diversity analysis,
#'     transformations preserve relative abundances while reducing dominance effects
#'   }
#'
#'   \item{\strong{Percentage/Coverage ("percent"):}}{
#'     Transformation: Arcsine-square root only
#'
#'     \strong{Rationale:} Data bounded between 0-1 (proportions) or 0-100 (percentages)
#'     have variance that depends on the mean (heteroscedastic). Arcsine-sqrt
#'     stabilizes variance across the range. This is the classical transformation
#'     for binomial proportions.
#'
#'     \strong{Input format:} Data must be proportions in the range of 0 to 1 OR percentages between 0 and 100.
#'     The function automatically detects and converts percentages to proportions.
#'
#'     \strong{Typical variables:} Percent cover, percent organic matter, percent
#'     sand/silt/clay, survival rates, colonization rates
#'
#'     \strong{Modern alternative:} Logit transformation or beta regression may be
#'     preferred for modern GLM approaches (Warton & Hui 2011), but arcsine-sqrt
#'     remains widely used and accepted
#'   }
#'
#'   \item{\strong{Presence-Absence ("pa"):}}{
#'     Transformation: Convert to binary (0/1)
#'
#'     \strong{Rationale:} Reduces data to occurrence only, removing abundance
#'     information. Useful when: (1) abundance data are unreliable, (2) analyzing
#'     distributional patterns, (3) calculating presence-based diversity indices
#'     (Jaccard, Sorensen), (4) reducing effects of sampling effort differences
#'
#'     \strong{Output:} All values > 0 become 1, all values = 0 remain 0.
#'     No normality testing performed as binary data are not continuous.
#'
#'     \strong{Use for:} Beta diversity (betapart package), incidence-based analyses,
#'     distribution modeling, occupancy analysis
#'   }
#' }
#'
#' \strong{Transformation Methods - Technical Details:}
#'
#' \describe{
#'   \item{\strong{Square root:} \eqn{\sqrt{x}}}{
#'     Stabilizes variance in count data (Poisson), moderate variance-stabilizing
#'     effect, suitable for low to moderate skewness
#'   }
#'
#'   \item{\strong{Log(x+1):} \eqn{\log(x + 1)}}{
#'     Strong variance stabilization, converts multiplicative to additive relationships,
#'     appropriate for lognormally distributed data, +1 handles zeros
#'   }
#'
#'   \item{\strong{Fourth root:} \eqn{x^{0.25}}}{
#'     Intermediate between sqrt and log, preferred in marine ecology for species
#'     abundance, maintains community structure while reducing dominance
#'   }
#'
#'   \item{\strong{Arcsine-sqrt:} \eqn{\arcsin(\sqrt{x})}}{
#'     Variance-stabilizing for binomial proportions, stretches ends between 0 to 1,
#'     compresses middle, formula: asin(sqrt(x)) where x ∈ (0,1)
#'   }
#'
#'   \item{\strong{Chord:} \eqn{x_{ij} / \sqrt{\sum x_{ij}^2}}}{
#'     Normalizes samples to unit length, preserves Euclidean distances,
#'     compositional standardization
#'   }
#'
#'   \item{\strong{Z-score:} \eqn{(x - \bar{x}) / \sigma}}{
#'     Centers to mean=0, scales to sd=1, makes variables comparable,
#'     does not change distribution shape
#'   }
#' }
#'
#' \strong{Output Usage:}
#' \itemize{
#'   \item \strong{normal_data:} Use for parametric tests (t-test, ANOVA,
#'     linear regression, PCA, RDA, discriminant analysis)
#'   \item \strong{non_normal_data:} Use for non-parametric tests (Mann-Whitney,
#'     Kruskal-Wallis, Spearman, NMDS, PERMANOVA, ANOSIM)
#' }
#'
#' @examples
#' # Example 1: Environmental variables
#' env_data <- data.frame(
#'   Station = paste0("St", 1:20),
#'   Temperature = rnorm(20, 15, 3),
#'   Depth = rexp(20, 0.1),
#'   Salinity = rnorm(20, 35, 2),
#'   Chlorophyll = rexp(20, 0.5)
#' )
#' result_env <- TransformData(env_data, data_type = "env")
#' result_env$transformation_report
#'
#' # Example 2: Biological abundance data
#' bio_data <- data.frame(
#'   Station = paste0("St", 1:20),
#'   Capitella = rpois(20, 50),
#'   Owenia = rpois(20, 10),
#'   Nephtys = rpois(20, 5)
#' )
#' result_bio <- TransformData(bio_data, data_type = "bio")
#'
#' # Example 3: Percentage data (0-100)
#' percent_data <- data.frame(
#'   Station = paste0("St", 1:20),
#'   Organic_Matter = runif(20, 0, 100),
#'   Sand = runif(20, 0, 100),
#'   Silt = runif(20, 0, 100)
#' )
#' result_percent <- TransformData(percent_data, data_type = "percent")
#'
#' # Example 4: Proportion data (0-1)
#' prop_data <- data.frame(
#'   Station = paste0("St", 1:20),
#'   Cover_Algae = runif(20, 0, 1),
#'   Cover_Coral = runif(20, 0, 1)
#' )
#' result_prop <- TransformData(prop_data, data_type = "percent")
#'
#' # Example 5: Convert to presence-absence
#' result_pa <- TransformData(bio_data, data_type = "pa")
#' head(result_pa$normal_data)  # All values are 0 or 1
#'
#' @references
#' Legendre, P. & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#'
#' Clarke, K.R. & Warwick, R.M. (2001). Change in Marine Communities: An Approach
#' to Statistical Analysis and Interpretation (2nd ed.). PRIMER-E, Plymouth.
#'
#' Field, J.G., Clarke, K.R. & Warwick, R.M. (1982). A practical strategy for
#' analysing multispecies distribution patterns. Marine Ecology Progress Series,
#' 8, 37-52. https://doi.org/10.3354/meps008037
#'
#' Zuur, A.F., Ieno, E.N. & Elphick, C.S. (2010). A protocol for data exploration
#' to avoid common statistical problems. Methods in Ecology and Evolution, 1(1), 3-14.
#' https://doi.org/10.1111/j.2041-210X.2009.00001.x
#'
#' Legendre, P. & Gallagher, E.D. (2001). Ecologically meaningful transformations
#' for ordination of species data. Oecologia, 129(2), 271-280.
#' https://doi.org/10.1007/s004420100716
#'
#' Warton, D.I. & Hui, F.K.C. (2011). The arcsine is asinine: the analysis of
#' proportions in ecology. Ecology, 92(1), 3-10. https://doi.org/10.1890/10-0340.1
#'
#' Sokal, R.R. & Rohlf, F.J. (1995). Biometry: The Principles and Practice of
#' Statistics in Biological Research (3rd ed.). W.H. Freeman.
#'
#' @seealso \code{\link{NormalTest}} for testing normality without transformation,
#'   \code{\link{PreData}} for data format preparation
#'
#' @importFrom stats shapiro.test ks.test sd
#' @export
TransformData <- function(DF, data_type = "env", alpha = 0.05, min_sd = 0, test = "shapiro") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!data_type %in% c("env", "bio", "percent", "pa")) {
    stop("data_type must be 'env', 'bio', 'percent', or 'pa'")
  }

  if (ncol(DF) < 2) {
    stop("DF must have at least 2 columns (names + variables)")
  }

  # Set row names from first column and remove it
  station_names <- as.character(DF[, 1])
  DF_numeric <- DF[, -1, drop = FALSE]
  rownames(DF_numeric) <- station_names

  # Select only numeric columns
  numeric_cols <- sapply(DF_numeric, is.numeric)
  if (sum(numeric_cols) == 0) {
    stop("DF must contain at least one numeric column")
  }
  DF_numeric <- DF_numeric[, numeric_cols, drop = FALSE]

  # Special case: Presence-Absence
  if (data_type == "pa") {
    pa_data <- apply_pa(DF_numeric)

    transformation_log <- data.frame(
      Variable = colnames(pa_data),
      Transformation = "presence-absence",
      Initial_p = NA,
      Final_p = NA,
      Status = "Converted to PA",
      stringsAsFactors = FALSE
    )

    summary_stats <- data.frame(
      Metric = c("Total variables", "Converted to presence-absence"),
      Count = c(ncol(pa_data), ncol(pa_data))
    )

    return(list(
      normal_data = pa_data,
      non_normal_data = data.frame(matrix(ncol = 0, nrow = nrow(pa_data))),
      transformation_report = transformation_log,
      summary = summary_stats,
      normal_vars = colnames(pa_data),
      non_normal_vars = character(0)
    ))
  }

  # Special case: Percentage/Proportion data
  if (data_type == "percent") {
    # Check if data are percentages between 0 and 100 or proportions between 0 and 1
    max_vals <- sapply(DF_numeric, max, na.rm = TRUE)

    # Convert percentages to proportions if needed
    DF_for_transform <- DF_numeric
    converted_vars <- character(0)

    for (var in names(max_vals)) {
      if (max_vals[var] > 1) {
        # Assume percentages, convert to proportions
        DF_for_transform[[var]] <- DF_numeric[[var]] / 100
        converted_vars <- c(converted_vars, var)
      }
    }

    if (length(converted_vars) > 0) {
      message("Note: The following variables were detected as percentages [0,100] and converted to proportions [0,1]: ",
              paste(converted_vars, collapse = ", "))
    }

    # Validate all values are in range 0 and 1
    if (any(DF_for_transform < 0 | DF_for_transform > 1, na.rm = TRUE)) {
      stop("All values must be proportions [0,1] or percentages [0,100]")
    }

    # Apply arcsine-sqrt transformation
    arcsin_data <- apply_arcsin(DF_for_transform)
    arcsin_results <- test_all_variables(arcsin_data, test)

    normal_vars <- arcsin_results$normal_vars
    non_normal_vars <- setdiff(colnames(arcsin_data), normal_vars)

    transformation_log <- data.frame(
      Variable = colnames(arcsin_data),
      Transformation = "arcsine-sqrt",
      Initial_p = NA,
      Final_p = arcsin_results$p_values,
      Status = ifelse(colnames(arcsin_data) %in% normal_vars, "Normal", "Non-normal"),
      stringsAsFactors = FALSE
    )

    transformation_log$Final_p <- round(transformation_log$Final_p, 4)

    normal_data <- if (length(normal_vars) > 0) {
      arcsin_data[, normal_vars, drop = FALSE]
    } else {
      data.frame(matrix(ncol = 0, nrow = nrow(arcsin_data)))
    }

    non_normal_data <- if (length(non_normal_vars) > 0) {
      DF_numeric[, non_normal_vars, drop = FALSE]
    } else {
      data.frame(matrix(ncol = 0, nrow = nrow(DF_numeric)))
    }

    summary_stats <- data.frame(
      Metric = c("Total variables", "Normalized by arcsine-sqrt", "Remaining non-normal"),
      Count = c(ncol(arcsin_data), length(normal_vars), length(non_normal_vars))
    )

    return(list(
      normal_data = normal_data,
      non_normal_data = non_normal_data,
      transformation_report = transformation_log,
      summary = summary_stats,
      normal_vars = normal_vars,
      non_normal_vars = non_normal_vars
    ))
  }

  # For "env" and "bio" data types - full transformation workflow

  # Step 1: Filter by standard deviation
  sds <- sapply(DF_numeric, sd, na.rm = TRUE)
  keep_vars <- names(sds)[sds > min_sd]

  if (length(keep_vars) == 0) {
    stop("No variables with SD > min_sd")
  }

  DF_filtered <- DF_numeric[, keep_vars, drop = FALSE]
  excluded_low_sd <- setdiff(names(DF_numeric), keep_vars)

  # Initialize tracking
  all_vars <- colnames(DF_filtered)
  original_data <- DF_filtered
  transformation_log <- data.frame(
    Variable = all_vars,
    Transformation = "none",
    Initial_p = NA,
    Final_p = NA,
    Status = "pending",
    stringsAsFactors = FALSE
  )

  # Containers for results
  normal_dfs <- list()
  current_non_normal <- DF_filtered

  # Step 2: Initial normality test
  initial_results <- test_all_variables(current_non_normal, test)

  # Separate normal and non-normal
  initially_normal <- initial_results$normal_vars
  if (length(initially_normal) > 0) {
    normal_dfs[["none"]] <- current_non_normal[, initially_normal, drop = FALSE]
    transformation_log$Transformation[transformation_log$Variable %in% initially_normal] <- "none"
    transformation_log$Initial_p[transformation_log$Variable %in% initially_normal] <-
      initial_results$p_values[initially_normal]
    transformation_log$Final_p[transformation_log$Variable %in% initially_normal] <-
      initial_results$p_values[initially_normal]
    transformation_log$Status[transformation_log$Variable %in% initially_normal] <- "Normal"

    current_non_normal <- current_non_normal[, setdiff(colnames(current_non_normal), initially_normal), drop = FALSE]
  }

  # Record initial p-values
  transformation_log$Initial_p[transformation_log$Variable %in% colnames(current_non_normal)] <-
    initial_results$p_values[colnames(current_non_normal)]

  # Define transformation sequence based on data type
  if (data_type == "env") {
    # Environmental: sqrt → log → fourth → chord → scale
    transform_sequence <- list(
      list(name = "sqrt", func = apply_sqrt),
      list(name = "log", func = apply_log),
      list(name = "fourth", func = apply_fourth),
      list(name = "chord", func = apply_chord),
      list(name = "z-score", func = apply_scale)
    )
  } else if (data_type == "bio") {
    # Biological: sqrt → fourth → log → chord → scale
    transform_sequence <- list(
      list(name = "sqrt", func = apply_sqrt),
      list(name = "fourth", func = apply_fourth),
      list(name = "log", func = apply_log),
      list(name = "chord", func = apply_chord),
      list(name = "z-score", func = apply_scale)
    )
  }

  # Apply transformations in sequence
  for (transform in transform_sequence) {
    if (ncol(current_non_normal) == 0) break

    transformed_data <- transform$func(current_non_normal)
    transform_results <- test_all_variables(transformed_data, test)

    normal_vars <- transform_results$normal_vars[transform_results$p_values[transform_results$normal_vars] > alpha]

    if (length(normal_vars) > 0) {
      normal_dfs[[transform$name]] <- transformed_data[, normal_vars, drop = FALSE]
      transformation_log$Transformation[transformation_log$Variable %in% normal_vars] <- transform$name
      transformation_log$Final_p[transformation_log$Variable %in% normal_vars] <-
        transform_results$p_values[normal_vars]
      transformation_log$Status[transformation_log$Variable %in% normal_vars] <- "Normal"

      current_non_normal <- current_non_normal[, setdiff(colnames(current_non_normal), normal_vars), drop = FALSE]
    }
  }

  # Mark remaining as non-normal
  if (ncol(current_non_normal) > 0) {
    remaining_vars <- colnames(current_non_normal)
    transformation_log$Status[transformation_log$Variable %in% remaining_vars] <- "Non-normal"
    final_test <- test_all_variables(current_non_normal, test)
    transformation_log$Final_p[transformation_log$Variable %in% remaining_vars] <-
      final_test$p_values[remaining_vars]
  }

  # Combine all normal data
  if (length(normal_dfs) > 0) {
    normal_data <- do.call(cbind, normal_dfs)
  } else {
    normal_data <- data.frame(matrix(ncol = 0, nrow = nrow(DF_filtered)))
    rownames(normal_data) <- rownames(DF_filtered)
  }

  # Non-normal data (original values)
  non_normal_vars <- transformation_log$Variable[transformation_log$Status == "Non-normal"]
  if (length(non_normal_vars) > 0) {
    non_normal_data <- original_data[, non_normal_vars, drop = FALSE]
  } else {
    non_normal_data <- data.frame(matrix(ncol = 0, nrow = nrow(DF_filtered)))
    rownames(non_normal_data) <- rownames(DF_filtered)
  }

  # Round p-values
  transformation_log$Initial_p <- round(transformation_log$Initial_p, 4)
  transformation_log$Final_p <- round(transformation_log$Final_p, 4)

  # Sort by status
  transformation_log <- transformation_log[order(transformation_log$Status, decreasing = TRUE), ]
  rownames(transformation_log) <- NULL

  # Create summary
  transform_counts <- table(transformation_log$Transformation[transformation_log$Status == "Normal"])

  summary_stats <- data.frame(
    Metric = c("Total variables", "Excluded (low SD)", "Initially normal",
               names(transform_counts),
               "Remaining non-normal", "Total normalized"),
    Count = c(
      length(all_vars) + length(excluded_low_sd),
      length(excluded_low_sd),
      sum(transformation_log$Transformation == "none" & transformation_log$Status == "Normal"),
      as.numeric(transform_counts),
      sum(transformation_log$Status == "Non-normal"),
      ncol(normal_data)
    )
  )

  # Return results
  return(list(
    normal_data = normal_data,
    non_normal_data = non_normal_data,
    transformation_report = transformation_log,
    summary = summary_stats,
    normal_vars = colnames(normal_data),
    non_normal_vars = colnames(non_normal_data)
  ))
}

# Helper functions (same as before)

test_all_variables <- function(DF, test) {
  p_values <- sapply(colnames(DF), function(var) {
    x <- DF[[var]]
    x_clean <- x[!is.na(x)]

    if (length(x_clean) < 3) return(NA)

    if (test == "shapiro") {
      if (length(x_clean) > 5000) return(NA)
      tryCatch(shapiro.test(x_clean)$p.value, error = function(e) NA)
    } else {
      mean_x <- mean(x_clean)
      sd_x <- sd(x_clean)
      tryCatch(ks.test(x_clean, "pnorm", mean = mean_x, sd = sd_x)$p.value,
               error = function(e) NA)
    }
  })

  normal_vars <- names(p_values)[!is.na(p_values) & p_values > 0.05]

  list(p_values = p_values, normal_vars = normal_vars)
}

apply_sqrt <- function(DF) {
  as.data.frame(lapply(DF, function(x) sqrt(x)))
}

apply_log <- function(DF) {
  as.data.frame(lapply(DF, function(x) log(x + 1)))
}

apply_fourth <- function(DF) {
  as.data.frame(lapply(DF, function(x) x^0.25))
}

apply_arcsin <- function(DF) {
  as.data.frame(lapply(DF, function(x) asin(sqrt(x))))
}

apply_pa <- function(DF) {
  as.data.frame(lapply(DF, function(x) as.numeric(x > 0)))
}

apply_chord <- function(DF) {
  row_norms <- sqrt(rowSums(DF^2, na.rm = TRUE))
  as.data.frame(sweep(DF, 1, row_norms, "/"))
}

apply_scale <- function(DF) {
  as.data.frame(scale(DF))
}


#' Prepare Data for Analysis
#'
#' Prepares community data by removing columns, setting row names, and optionally
#' transposing the data frame. Useful for cleaning and reformatting data before
#' diversity analysis.
#'
#' @param DF A data.frame containing community data.
#' @param ReCols Integer. Number of columns to remove from the beginning of the
#'   data frame (default = 0). Useful when the first column(s) contain metadata
#'   that should be excluded.
#' @param NCols Integer. Column number to use as row names (default = 1).
#'   This column will be converted to row names and then removed from the data.
#' @param transpose Logical. Should the data frame be transposed after processing?
#'   (default = FALSE). Set to TRUE to switch rows and columns (e.g., convert
#'   species-in-rows to species-in-columns format).
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item Row names set from the specified column
#'     \item Specified columns removed
#'     \item Optionally transposed
#'   }
#'
#' @details
#' This function streamlines common data preparation tasks:
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item Remove first \code{ReCols} columns (if ReCols > 0)
#'   \item Set row names from column \code{NCols}
#'   \item Remove the names column
#'   \item Transpose if \code{transpose = TRUE}
#' }
#'
#' \strong{Common use cases:}
#' \itemize{
#'   \item Converting Excel/CSV imports to analysis-ready format
#'   \item Handling data with multiple ID columns
#'   \item Switching between wide and long formats
#'   \item Preparing data for diversity functions
#' }
#'
#' @examples
#' # Example 1: Simple case - set row names from column 1
#' data_raw <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 12, 3),
#'   Owenia = c(23, 18, 45),
#'   Nephtys = c(12, 34, 23)
#' )
#' PreData(data_raw, NCols = 1, transpose = FALSE)
#'
#' # Example 2: Remove first column, use second as names, and transpose
#' data_raw2 <- data.frame(
#'   ID = 1:3,
#'   Station = c("St1", "St2", "St3"),
#'   Capitella = c(45, 12, 3),
#'   Owenia = c(23, 18, 45)
#' )
#' PreData(data_raw2, ReCols = 1, NCols = 1, transpose = TRUE)
#'
#' # Example 3: Transpose without removing columns
#' data_raw3 <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys"),
#'   St1 = c(45, 23, 12),
#'   St2 = c(12, 18, 34),
#'   St3 = c(3, 45, 23)
#' )
#' PreData(data_raw3, NCols = 1, transpose = TRUE)
#'
#' @importFrom stats shapiro.test ks.test sd
#' @export
PreData <- function(DF, ReCols = 0, NCols = 1, transpose = FALSE) {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!is.numeric(ReCols) || ReCols < 0 || ReCols >= ncol(DF)) {
    stop("ReCols must be a non-negative integer less than the number of columns")
  }

  if (!is.numeric(NCols) || NCols < 1 || NCols > ncol(DF)) {
    stop("NCols must be a positive integer within the number of columns")
  }

  if (!is.logical(transpose)) {
    stop("transpose must be TRUE or FALSE")
  }

  # Remove initial columns if specified
  if (ReCols > 0) {
    DF <- DF[, -(1:ReCols), drop = FALSE]
  }

  # Set row names from specified column
  rownames(DF) <- DF[, NCols]

  # Remove the names column
  DF <- DF[, -NCols, drop = FALSE]

  # Transpose if requested
  if (transpose) {
    DF <- as.data.frame(t(DF))
  }

  return(DF)
}


#' Test Multiple Variables for Normality
#'
#' Performs normality tests (Shapiro-Wilk and/or Kolmogorov-Smirnov) on all numeric
#' columns of a data frame and returns a formatted summary table with test results
#' and recommendations for transformation.
#'
#' @param DF A data.frame containing numeric variables to test for normality.
#'   The first column should contain station names (if shape = "w") or species
#'   names (if shape = "l").
#' @param tests Character vector specifying which test(s) to perform:
#'   \itemize{
#'     \item \code{"shapiro"} - Shapiro-Wilk test only (default)
#'     \item \code{"ks"} - Kolmogorov-Smirnov test only
#'     \item \code{"both"} - Both tests
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns
#'       (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns
#'       (first column = species names)
#'   }
#'
#' @return A data.frame with results depending on the test(s) selected.
#'   See Details for column descriptions.
#'
#' @details
#' \strong{Shapiro-Wilk Test:}
#'
#' The Shapiro-Wilk test is one of the most powerful normality tests, particularly
#' effective for small to moderate sample sizes (3 <= n <= 5000).
#'
#' \strong{Kolmogorov-Smirnov Test:}
#'
#' The Kolmogorov-Smirnov test compares the empirical cumulative distribution
#' function to the theoretical normal distribution.
#'
#' \strong{Choosing Between Tests:}
#' \itemize{
#'   \item \strong{Small samples (n < 50):} Use Shapiro-Wilk (more powerful)
#'   \item \strong{Large samples (n > 5000):} Use Kolmogorov-Smirnov
#'   \item \strong{Moderate samples (50 <= n <= 5000):} Use both for confirmation
#' }
#'
#' @examples
#' # Example 1: Test with Shapiro-Wilk (default)
#' data_wide <- data.frame(
#'   Station = c("St1", "St2", "St3"),
#'   Var1 = rnorm(3),
#'   Var2 = rexp(3)
#' )
#' NormalTest(data_wide, tests = "shapiro", shape = "w")
#'
#' # Example 2: Long format
#' data_long <- data.frame(
#'   Species = c("Sp1", "Sp2"),
#'   St1 = rnorm(2),
#'   St2 = rnorm(2)
#' )
#' NormalTest(data_long, tests = "shapiro", shape = "l")

#'
#' @references
#' Shapiro, S.S. & Wilk, M.B. (1965). An analysis of variance test for normality
#' (complete samples). Biometrika, 52(3-4), 591-611.
#' https://doi.org/10.1093/biomet/52.3-4.591
#'
#' Royston, P. (1982). An extension of Shapiro and Wilk's W test for normality
#' to large samples. Applied Statistics, 31(2), 115-124.
#'
#' Kolmogorov, A. (1933). Sulla determinazione empirica di una legge di
#' distribuzione. Giornale dell'Istituto Italiano degli Attuari, 4, 83-91.
#'
#' Razali, N.M. & Wah, Y.B. (2011). Power comparisons of Shapiro-Wilk,
#' Kolmogorov-Smirnov, Lilliefors and Anderson-Darling tests. Journal of
#' Statistical Modeling and Analytics, 2(1), 21-33.
#'
#' @seealso \code{\link[stats]{shapiro.test}} for Shapiro-Wilk test,
#'   \code{\link[stats]{ks.test}} for Kolmogorov-Smirnov test
#' @importFrom stats shapiro.test ks.test sd
#' @export
NormalTest <- function(DF, tests = "shapiro", shape = "w") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!tests %in% c("shapiro", "ks", "both")) {
    stop("tests must be 'shapiro', 'ks', or 'both'")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    Variables <- DF[, 1]
    DF_numeric <- as.data.frame(t(DF[, -1]))
    colnames(DF_numeric) <- Variables
    rownames(DF_numeric) <- colnames(DF[, -1])
  } else {
    # Wide format: Stations in first column, variables in other columns
    DF_numeric <- DF[, -1]
    Variables <- colnames(DF_numeric)
    rownames(DF_numeric) <- DF[, 1]
  }

  # Select only numeric columns
  numeric_cols <- sapply(DF_numeric, is.numeric)

  if (sum(numeric_cols) == 0) {
    stop("DF must contain at least one numeric column")
  }

  DF_numeric <- DF_numeric[, numeric_cols, drop = FALSE]

  # Initialize results list
  results_list <- list()

  # Shapiro-Wilk test
  if (tests %in% c("shapiro", "both")) {

    shapiro_results <- lapply(DF_numeric, function(x) {
      x_clean <- x[!is.na(x)]

      if (length(x_clean) < 3) {
        return(list(statistic = NA, p.value = NA, note = "Insufficient data (n < 3)"))
      }

      if (length(x_clean) > 5000) {
        return(list(statistic = NA, p.value = NA, note = "Sample too large (n > 5000)"))
      }

      tryCatch(
        shapiro.test(x_clean),
        error = function(e) list(statistic = NA, p.value = NA, note = as.character(e))
      )
    })

    W_stats <- sapply(shapiro_results, function(x) {
      if (is.list(x) && "statistic" %in% names(x)) {
        return(as.numeric(x$statistic))
      } else {
        return(NA)
      }
    })

    p_SW <- sapply(shapiro_results, function(x) {
      if (is.list(x) && "p.value" %in% names(x)) {
        return(as.numeric(x$p.value))
      } else {
        return(NA)
      }
    })

    results_list$W <- round(W_stats, 3)
    results_list$p_SW <- round(p_SW, 3)
    results_list$Status_SW <- ifelse(
      is.na(p_SW),
      "Unable to test",
      ifelse(p_SW > 0.05, "Normal", "Transform")
    )
  }

  # Kolmogorov-Smirnov test
  if (tests %in% c("ks", "both")) {

    ks_results <- lapply(DF_numeric, function(x) {
      x_clean <- x[!is.na(x)]

      if (length(x_clean) < 3) {
        return(list(statistic = NA, p.value = NA))
      }

      # Calculate mean and SD for comparison
      mean_x <- mean(x_clean)
      sd_x <- sd(x_clean)

      tryCatch(
        ks.test(x_clean, "pnorm", mean = mean_x, sd = sd_x),
        error = function(e) list(statistic = NA, p.value = NA)
      )
    })

    D_stats <- sapply(ks_results, function(x) {
      if (is.list(x) && "statistic" %in% names(x)) {
        return(as.numeric(x$statistic))
      } else {
        return(NA)
      }
    })

    p_KS <- sapply(ks_results, function(x) {
      if (is.list(x) && "p.value" %in% names(x)) {
        return(as.numeric(x$p.value))
      } else {
        return(NA)
      }
    })

    results_list$D <- round(D_stats, 3)
    results_list$p_KS <- round(p_KS, 3)
    results_list$Status_KS <- ifelse(
      is.na(p_KS),
      "Unable to test",
      ifelse(p_KS > 0.05, "Normal", "Transform")
    )
  }

  # Create output data frame
  output <- data.frame(
    Variables = names(DF_numeric),
    stringsAsFactors = FALSE
  )

  # Add results based on tests performed
  if (tests == "shapiro") {
    output$W <- results_list$W
    output$p_SW <- results_list$p_SW
    output$Status <- results_list$Status_SW

  } else if (tests == "ks") {
    output$D <- results_list$D
    output$p_KS <- results_list$p_KS
    output$Status <- results_list$Status_KS

  } else {  # both
    output$W <- results_list$W
    output$p_SW <- results_list$p_SW
    output$D <- results_list$D
    output$p_KS <- results_list$p_KS

    # Combined recommendation
    output$Recommendation <- ifelse(
      results_list$Status_SW == "Normal" & results_list$Status_KS == "Normal",
      "Normal (both tests)",
      ifelse(
        results_list$Status_SW == "Transform" & results_list$Status_KS == "Transform",
        "Transform (both tests agree)",
        ifelse(
          results_list$Status_SW == "Unable to test" | results_list$Status_KS == "Unable to test",
          "Check data quality",
          "Inconclusive (tests disagree - visual inspection recommended)"
        )
      )
    )
  }

  # Sort by status/recommendation
  if ("Recommendation" %in% colnames(output)) {
    output <- output[order(output$Recommendation), ]
  } else {
    output <- output[order(output$Status), ]
  }

  rownames(output) <- NULL

  return(output)
}


#' Filter Rare Species from Community Data
#'
#' Removes rare species from community data based on occurrence frequency
#' and/or total abundance thresholds. Useful for reducing noise and focusing
#' analyses on more common, ecologically relevant species.
#'
#' @param DF A data.frame containing community data. The first column should
#'   contain station names (if shape = "w") or species names (if shape = "l").
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, species in columns
#'       (first column = station names) (default)
#'     \item \code{"l"} - Long format: species in rows, stations in columns
#'       (first column = species names)
#'   }
#' @param mocc  Minimum proportion of samples where species must occur
#'   (default = 0.05, i.e., 5\% of samples). Species occurring in fewer samples
#'   are removed. Set to 0 to disable occurrence filtering.
#' @param mab Minimum total abundance across all samples (default = 1).
#'   Species with total abundance below this threshold are removed. Set to 0
#'   to disable abundance filtering.
#' @param method Filtering method (default = "and"):
#'   \itemize{
#'     \item \code{"and"} - Species must meet BOTH criteria (occurrence AND abundance)
#'     \item \code{"or"} - Species must meet EITHER criterion (occurrence OR abundance)
#'   }
#'
#' @return A list with four elements:
#'   \itemize{
#'     \item \code{filtered_data} - Data frame with rare species removed, same
#'       format as input (stations as rows if shape="w")
#'     \item \code{removed_species} - Character vector of removed species names
#'     \item \code{kept_species} - Character vector of retained species names
#'     \item \code{summary} - Data frame with filtering statistics including
#'       original species count, removed count, retained count, and criteria used
#'   }
#'
#' @details
#' Rare species are common in ecological datasets but can:
#' \itemize{
#'   \item Add noise to multivariate analyses (ordination, clustering)
#'   \item Inflate apparent beta diversity
#'   \item Be artifacts of sampling error or misidentification
#'   \item Have disproportionate influence in some diversity indices
#' }
#'
#' \strong{Filtering Criteria:}
#'
#' \describe{
#'   \item{\strong{Occurrence (mocc):}}{
#'     Proportion of samples where species is present (abundance > 0).
#'
#'     Example: mocc = 0.05 with 100 samples means species must
#'     occur in at least 5 samples.
#'
#'     \strong{Use when:} Focusing on widespread species, reducing singleton effects,
#'     preparing data for presence-absence analyses
#'   }
#'
#'   \item{\strong{Abundance (mab):}}{
#'     Total number of individuals across all samples.
#'
#'     Example: mab = 10 means species must have at least 10 total
#'     individuals across all samples combined.
#'
#'     \strong{Use when:} Removing species with very few individuals that may be
#'     sampling artifacts, focusing on numerically important species
#'   }
#'
#'   \item{\strong{Method "and" (default):}}{
#'     Species must satisfy BOTH occurrence AND abundance thresholds.
#'     More conservative - removes more species.
#'
#'     \strong{Use when:} You want only well-represented species
#'   }
#'
#'   \item{\strong{Method "or":}}{
#'     Species must satisfy EITHER occurrence OR abundance threshold.
#'     More permissive - removes fewer species.
#'
#'     \strong{Use when:} You want to keep species that are either widespread
#'     (even if not abundant) OR abundant (even if localized)
#'   }
#' }
#'
#' \strong{Recommended Thresholds:}
#' \itemize{
#'   \item \strong{Exploratory analysis:} mocc = 0.05 (5\%),
#'     mab = 5-10
#'   \item \strong{Ordination (NMDS, PCA):} mocc = 0.10 (10\%),
#'     mab = 10
#'   \item \strong{Indicator species:} mocc = 0.20 (20\%),
#'     mab = 20
#'   \item \strong{Diversity indices:} Use with caution - may bias richness
#' }
#'
#' \strong{Important Considerations:}
#' \itemize{
#'   \item Filtering affects richness-based diversity indices (reduce values)
#'   \item May remove ecologically important but rare species (endangered, specialists)
#'   \item Should be justified and reported in methods
#'   \item Consider ecological context - a "rare" predator may be ecologically important
#'   \item For beta diversity, filtering can reduce artifactual dissimilarity
#' }
#'
#' @examples
#' # Example 1: Basic filtering - wide format
#' community_data <- data.frame(
#'   Station = paste0("St", 1:20),
#'   Common_sp1 = rpois(20, 50),      # Common: occurs in all stations
#'   Rare_sp1 = c(rep(0, 18), 3, 2),  # Rare: only 2 stations, low abundance
#'   Common_sp2 = rpois(20, 30),      # Common
#'   Rare_sp2 = c(1, rep(0, 19))      # Very rare: only 1 station
#' )
#'
#' # Filter: must occur in 10% of samples AND have 10+ total individuals
#' result <- FilterRare(community_data, shape = "w",
#'                      mocc = 0.10, mab = 10)
#' result$summary
#' result$removed_species
#'
#' # Example 2: Long format
#' community_long <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys", "Rare_sp"),
#'   St1 = c(45, 23, 12, 0),
#'   St2 = c(50, 18, 15, 1),
#'   St3 = c(48, 25, 10, 0)
#' )
#' result2 <- FilterRare(community_long, shape = "l",
#'                       mocc = 0.30, mab = 5)
#'
#' # Example 3: Only occurrence-based filtering
#' result3 <- FilterRare(community_data, shape = "w",
#'                       mocc = 0.15, mab = 0)
#'
#' # Example 4: Only abundance-based filtering
#' result4 <- FilterRare(community_data, shape = "w",
#'                       mocc = 0, mab = 20)
#'
#' # Example 5: Permissive filtering (OR method)
#' result5 <- FilterRare(community_data, shape = "w",
#'                       mocc = 0.25, mab = 50,
#'                       method = "or")
#'
#' @references
#' Cao, Y., Williams, D.D. & Larsen, D.P. (2002). Comparison of ecological
#' communities: the problem of sample representativeness. Ecological Monographs,
#' 72(1), 41-56. https://doi.org/10.1890/0012-9615(2002)072 (0041:COECTP)2.0.CO;2
#'
#' Poos, M.S. & Jackson, D.A. (2012). Addressing the removal of rare species in
#' multivariate bioassessments: the impact of methodological choices. Ecological
#' Indicators, 18, 82-90. https://doi.org/10.1016/j.ecolind.2011.10.008
#'
#' Marchant, R. (2002). Do rare species have any place in multivariate analysis
#' for bioassessment? Journal of the North American Benthological Society, 21(2),
#' 311-313.
#'
#' Legendre, P. & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#'
#' @seealso \code{\link{PreData}} for data format preparation
#'
#' @export
FilterRare <- function(DF, shape = "w", mocc = 0.05,
                       mab = 1, method = "and") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'")
  }

  if (mocc < 0 || mocc > 1) {
    stop("mocc must be between 0 and 1")
  }

  if (mab < 0) {
    stop("mab must be non-negative")
  }

  if (!method %in% c("and", "or")) {
    stop("method must be 'and' or 'or'")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species in first column, stations in other columns
    species_names <- as.character(DF[, 1])
    abundance_matrix <- as.matrix(DF[, -1])
    rownames(abundance_matrix) <- species_names
    station_names <- colnames(DF[, -1])
  } else {
    # Wide format: Stations in first column, species in other columns
    station_names <- as.character(DF[, 1])
    abundance_matrix <- as.matrix(DF[, -1])
    rownames(abundance_matrix) <- station_names
    species_names <- colnames(abundance_matrix)
    # Transpose so species are in rows for filtering
    abundance_matrix <- t(abundance_matrix)
  }

  # Calculate filtering metrics for each species
  n_stations <- ncol(abundance_matrix)

  # Occurrence: proportion of samples where species is present
  occurrence <- rowSums(abundance_matrix > 0) / n_stations

  # Total abundance: sum across all samples
  total_abundance <- rowSums(abundance_matrix)

  # Determine which species to keep based on method
  if (method == "and") {
    # Both criteria must be met
    keep_species <- (occurrence >= mocc) & (total_abundance >= mab)
  } else {  # method == "or"
    # Either criterion can be met
    keep_species <- (occurrence >= mocc) | (total_abundance >= mab)
  }

  # Identify kept and removed species
  kept_species_names <- species_names[keep_species]
  removed_species_names <- species_names[!keep_species]

  # Filter the abundance matrix
  filtered_matrix <- abundance_matrix[keep_species, , drop = FALSE]

  # Convert back to original format
  if (shape == "l") {
    # Return as long format
    filtered_data <- data.frame(
      Species = rownames(filtered_matrix),
      filtered_matrix,
      row.names = NULL,
      check.names = FALSE
    )
    colnames(filtered_data) <- c(colnames(DF)[1], station_names)
  } else {
    # Return as wide format (transpose back)
    filtered_matrix <- t(filtered_matrix)
    filtered_data <- data.frame(
      Station = station_names,
      filtered_matrix,
      row.names = NULL,
      check.names = FALSE
    )
    colnames(filtered_data) <- c(colnames(DF)[1], kept_species_names)
  }

  # Create summary statistics
  summary_stats <- data.frame(
    Metric = c(
      "Original species count",
      "Species removed",
      "Species retained",
      "Removal percentage",
      "Min occurrence threshold",
      "Min abundance threshold",
      "Filter method"
    ),
    Value = c(
      length(species_names),
      length(removed_species_names),
      length(kept_species_names),
      round(100 * length(removed_species_names) / length(species_names), 1),
      mocc,
      mab,
      method
    ),
    stringsAsFactors = FALSE
  )

  # Return results
  return(list(
    filtered_data = filtered_data,
    removed_species = removed_species_names,
    kept_species = kept_species_names,
    summary = summary_stats
  ))
}


