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
#'   Special case: If NCols = 0, the function will automatically transpose the
#'   table, converting original column names into the first column (character)
#'   and the rest as numeric columns.
#' @param transpose Logical. Should the data frame be transposed after processing?
#'   (default = FALSE). Only applies when NCols > 0. Ignored when NCols = 0.
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item Row names set from the specified column (if NCols > 0)
#'     \item Specified columns removed
#'     \item Optionally transposed (if NCols > 0 and transpose = TRUE)
#'     \item If NCols = 0: automatically transposed data where original column
#'           names become the first column (character) and the rest as numeric columns
#'   }
#'
#' @details
#' This function streamlines common data preparation tasks:
#'
#' \strong{Workflow (NCols > 0 - DEFAULT):}
#' \enumerate{
#'   \item Remove first \code{ReCols} columns (if ReCols > 0)
#'   \item Set row names from column \code{NCols}
#'   \item Remove the names column
#'   \item Transpose if \code{transpose = TRUE}
#' }
#'
#' \strong{Special workflow (NCols = 0):}
#' \enumerate{
#'   \item Remove first \code{ReCols} columns (if ReCols > 0)
#'   \item Automatically transpose the data frame excluding the first column
#'   \item Set column names using the first column values (species names)
#'   \item Add a new first column with the original column names (sample IDs)
#'   \item Convert the new first column to character
#'   \item Convert all other columns to numeric
#' }
#' Note: When NCols = 0, transposition is automatic, the \code{transpose}
#' argument is ignored.
#'
#' \strong{Common use cases:}
#' \itemize{
#'   \item Converting Excel/CSV imports to analysis-ready format
#'   \item Handling data with multiple ID columns
#'   \item Switching between wide and long formats
#'   \item Preparing data for diversity functions
#'   \item Converting species-in-columns format with proper data types
#' }
#'
#' @examples
#' # Example 1: Simple case - set row names from column 1 (DEFAULT BEHAVIOR)
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
#' # Example 3: Transpose without removing columns (ORIGINAL BEHAVIOR)
#' data_raw3 <- data.frame(
#'   Species = c("Capitella", "Owenia", "Nephtys"),
#'   St1 = c(45, 23, 12),
#'   St2 = c(12, 18, 34),
#'   St3 = c(3, 45, 23)
#' )
#' PreData(data_raw3, NCols = 1, transpose = TRUE)
#'
#' # Example 4: NEW FEATURE - NCols = 0 auto-transposes with column names as first column
#'
#' @importFrom stats shapiro.test ks.test sd
#' @export
PreData <- function(DF, ReCols = 0, NCols = 1, transpose = FALSE) {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  # NEW FEATURE: Special case - NCols = 0
  # Logic equivalent to:
  # b = data.frame(t(Bio[,-1]))
  # colnames(b) = Bio$Species
  # b$var = rownames(b)
  # b = b[,c(dim(b)[2],1:dim(b)[2]-1)]
  if (NCols == 0) {
    # Remove initial columns if specified
    if (ReCols > 0) {
      if (ReCols >= ncol(DF)) {
        stop("ReCols must be less than the number of columns")
      }
      DF <- DF[, -(1:ReCols), drop = FALSE]
    }

    # Step 1: Transpose the data excluding the first column
    DF_transposed <- as.data.frame(t(DF[, -1, drop = FALSE]), stringsAsFactors = FALSE)

    # Step 2: Set column names using the first column of original data
    colnames(DF_transposed) <- DF[, 1]

    # Step 3: Add a new column with the original column names (sample IDs)
    DF_transposed$var <- rownames(DF_transposed)

    # Step 4: Reorder columns to put 'var' as the first column
    DF_transposed <- DF_transposed[, c(dim(DF_transposed)[2], 1:(dim(DF_transposed)[2] - 1))]

    # Step 5: Convert 'var' column to character
    DF_transposed[, 1] <- as.character(DF_transposed[, 1])

    # Step 6: Convert all other columns to numeric
    if (ncol(DF_transposed) > 1) {
      for (i in 2:ncol(DF_transposed)) {
        DF_transposed[, i] <- as.numeric(as.character(DF_transposed[, i]))
      }
    }

    return(DF_transposed)
  }

  # ORIGINAL BEHAVIOR: Normal case - NCols > 0 (DEFAULT)
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

#' Transform Data for Statistics or Multivariate Analysis
#'
#' @description
#' Transforms environmental and biological data to achieve normality.
#' Supports multiple data types with automatic or user-specified
#' transformations, ideal for preparing data before statistical analyses like
#' PCA, RDA, ANOVA.
#'
#' @param DF A data.frame containing the data to transform. Must include at least
#'   one identifier column and one or more numeric columns.
#' @param shape Character. Format of the input data: \code{"w"} (wide) or \code{"l"} (long).
#'   \itemize{
#'     \item \code{"w"}: Samples in first column, variables in subsequent columns (default)
#'     \item \code{"l"}: Variables in first column, samples in subsequent columns
#'   }
#' @param data_type Character. Type of data to transform. Options:
#'   \itemize{
#'     \item \code{"env"}: Environmental data - sequential transformations until normality
#'     \item \code{"bio"}: Biological/abundance data - single transformation
#'     \item \code{"percent"}: Percentage/proportion data - arcsine-square root
#'     \item \code{"pa"}: Presence-absence conversion
#'     \item \code{"scale"}: Z-score standardization only (no transformation)
#'   }
#'
#' @param version Character. Output format: \code{"l"} (long) or \code{"s"} (short).
#'   \itemize{
#'     \item \code{"l"}: Returns complete data.frame with first column preserved (default)
#'     \item \code{"s"}: Returns list with separated transformed/untransformed data and reports
#'   }
#' @param alpha Numeric. Significance level for normality tests (default = 0.05).
#'   Only used when \code{data_type = "env"}.
#' @param min_sd Numeric. Minimum standard deviation threshold. Variables with SD
#'   <= \code{min_sd} are excluded (default = 0). Only used when \code{data_type = "env"}.
#' @param test Character. Normality test to use: \code{"shapiro"} (Shapiro-Wilk) or
#'   \code{"ks"} (Kolmogorov-Smirnov). Default = \code{"shapiro"}.
#' @param trans Character. Transformation for biological data. Options:
#'   \itemize{
#'     \item \code{"r4"}: Fourth root (x^0.25) - default
#'     \item \code{"l1"}: Log(x + 1)
#'     \item \code{"sq"}: Square root
#'     \item \code{"Ch"}: Chord transformation
#'   }
#' @param scale Logical. If \code{TRUE}, applies z-score standardization to the
#'   final data. For \code{data_type = "env"}, scales all variables (normal + non-normal).
#'   Default = \code{FALSE}.
#'
#' @param rel_type Character. Type of relativization when data_type = "rel".
#'   Options: \code{"total"} (divide by row sum) or \code{"max"} (divide by row maximum).
#'   Default = \code{"total"}.
#'
#' @return
#' If \code{version = "l"}: A data.frame with the first column preserved and all
#' variables (transformed/untransformed). Transformation metadata stored as attributes.
#'
#' If \code{version = "s"}: A list containing:
#'   \itemize{
#'     \item \code{normal_data} / \code{transformed_data}: Transformed variables
#'     \item \code{non_normal_data}: Untransformed variables (only for \code{"env"})
#'     \item \code{transformation_report}: Detailed transformation log
#'     \item \code{summary}: Summary statistics
#'   }
#'
#' @details
#' \strong{Environmental data (\code{"env"}):}
#' Applies sequential transformations to each variable until normality is achieved
#' (p > alpha) or all options are exhausted:
#' \enumerate{
#'   \item No transformation (test original data)
#'   \item Square root
#'   \item Log(x + 1)
#'   \item Fourth root
#' }
#' Variables that remain non-normal are kept in their original scale.
#'
#' \strong{Biological data (\code{"bio"}):}
#' Applies a single transformation to all variables. Common choices:
#' \itemize{
#'   \item Fourth root: Recommended for count/abundance data
#'   \item Log(x+1): For highly skewed data with zeros
#'   \item Chord: For community composition data
#' }
#'
#' \strong{Percentage data (\code{"percent"}):}
#' Automatically converts percentages (>1) to proportions (0-1) and applies
#' arcsine-square root transformation. Useful for granulometry, cover, or
#' relative abundance data.
#'
#' \strong{Presence-Absence (\code{"pa"}):}
#' Converts all values > 0 to 1 and 0 to 0. Useful for community composition
#' analyses (Jaccard, Sørensen indices).
#'
#' \strong{Scale only (\code{"scale"}):}
#' Applies z-score standardization without any transformation. Centers data to
#' mean = 0 and scales to SD = 1.
#'
#' @note
#' When \code{version = "l"}, transformation metadata can be accessed via:
#' \code{attr(result, "summary")} and \code{attr(result, "transformation_report")}
#'
#' @seealso
#' \code{\link[vegan]{decostand}} for alternative transformations
#'
#' @examples
#'
#' # Example 1: Environmental data (sequential transformations)
#'
#' # Create environmental data with modified station codes and variance
#' set.seed(123)
#' env_data <- data.frame(
#'   Sample = c("S01", "S02", "S03", "S04", "S05", "S06", "S07",
#'              "M01", "M02", "M03", "M04", "M05",
#'              "P01", "P02", "P03", "P04", "P05"),
#'   Temperature = rnorm(17, mean = 18, sd = 3.5),
#'   Salinity = rnorm(17, mean = 35, sd = 2.8),
#'   DO = rnorm(17, mean = 6.5, sd = 1.2)
#' )
#'
#' # Apply environmental transformations
#' env_trans <- TransformData(env_data,
#'                            shape = "w",
#'                            data_type = "env",
#'                            version = "l")
#' head(env_trans)
#'
#' # View transformation summary
#' attr(env_trans, "summary")
#' attr(env_trans, "transformation_report")
#'
#' # With scaling and different alpha
#' env_scaled <- TransformData(env_data,
#'                             data_type = "env",
#'                             alpha = 0.01,
#'                             min_sd = 0.5,
#'                             scale = TRUE)
#'
#'
#' # Example 2: Biological data (single transformation)
#'
#' # Create biological abundance data with modified species names
#' set.seed(456)
#' bio_data <- data.frame(
#'   Station = c("A1", "A2", "A3", "B1", "B2", "B3", "B4", "B5",
#'               "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4", "D5"),
#'   sp1 = round(rlnorm(17, meanlog = 2, sdlog = 1.5)),
#'   sp2 = round(rlnorm(17, meanlog = 1.5, sdlog = 1.2)),
#'   sp3 = round(rlnorm(17, meanlog = 3, sdlog = 1.8)),
#'   sp4 = round(rlnorm(17, meanlog = 2.5, sdlog = 1.3)),
#'   sp5 = round(rlnorm(17, meanlog = 1, sdlog = 1.0))
#' )
#'
#' # Default: Fourth root transformation
#' bio_r4 <- TransformData(bio_data,
#'                         shape = "w",
#'                         data_type = "bio",
#'                         trans = "r4")
#' head(bio_r4)
#'
#' # Log(x+1) transformation
#' bio_log <- TransformData(bio_data,
#'                          data_type = "bio",
#'                          trans = "l1")
#'
#' # Chord transformation with scaling
#' bio_chord <- TransformData(bio_data,
#'                            data_type = "bio",
#'                            trans = "Ch",
#'                            scale = TRUE)
#'
#' # Square root transformation
#' bio_sqrt <- TransformData(bio_data,
#'                           data_type = "bio",
#'                           trans = "sq")
#'
#'
#' # Example 3: Long  format (species in first column)
#'
#' # Same data in long format
#' bio_long <- data.frame(
#'   Species = c("sp1", "sp2", "sp3", "sp4", "sp5"),
#'   A1 = c(12, 5, 45, 8, 23),
#'   A2 = c(8, 3, 67, 12, 15),
#'   B1 = c(34, 18, 23, 6, 42),
#'   B2 = c(56, 7, 89, 15, 8),
#'   C1 = c(23, 12, 34, 9, 17)
#' )
#'
#' bio_long_trans <- TransformData(bio_long,
#'                                 shape = "l",
#'                                 data_type = "bio",
#'                                 trans = "r4")
#' head(bio_long_trans)
#'
#' #' # Example 4: Percentage/Proportion data
#'
#' # Create percentage data (granulometry)
#' set.seed(789)
#' percent_data <- data.frame(
#' Site = c("X1", "X2", "X3", "Y1", "Y2", "Y3", "Z1", "Z2"),
#' Sand = round(runif(8, 40, 80), 1),
#' Silt = round(runif(8, 10, 35), 1),
#' Clay = round(runif(8, 5, 25), 1)
#' )
#'
#' # Ensure sum is approximately 100
#' percent_data2 <- data.frame(Site = percent_data$Site,
#'                             (percent_data[,-1] / rowSums(percent_data[,-1]) * 100))
#'
#' # Apply arcsine-square root transformation
#' percent_trans <- TransformData(percent_data2,
#'                                shape = "w",
#'                                data_type = "percent")
#' head(percent_trans)
#'
#' # With scaling
#' percent_scaled <- TransformData(percent_data2,
#'                                 data_type = "percent",
#'                                 scale = TRUE)
#'
#' # Example 5: Relativization (USA bio_data del Ejemplo 2)
#' bio_rel <- TransformData(bio_data,    # ← bio_data ya existe
#'                          data_type = "rel",
#'                          rel_type = "total")
#'
#' # Example 6: Presence-Absence conversion
#'
#' # Convert abundance to presence-absence
#' pa_data <- TransformData(bio_data,
#'                          data_type = "pa")
#' head(pa_data)
#'
#' # Example 7: Scale only (z-score standardization)
#'
#' # Apply only scaling without transformation
#' scaled_only <- TransformData(env_data,
#'                              data_type = "scale")
#' head(scaled_only)
#'
#' # Check: mean ~ 0, sd ~ 1
#' colMeans(scaled_only[, -1])
#' apply(scaled_only[, -1], 2, sd)
#'
#' # Example 8: Short version output
#'
#' # Get list output for more detailed access
#' env_list <- TransformData(env_data,
#'                           data_type = "env",
#'                           version = "s")
#'
#' env_list$summary
#' head(env_list$normal_data)
#' head(env_list$non_normal_data)
#' env_list$normal_vars
#' env_list$non_normal_vars
#'
#' # Example 9: Different normality tests
#'
#' # Using Kolmogorov-Smirnov test (less conservative)
#' env_ks <- TransformData(env_data,
#'                         data_type = "env",
#'                         test = "ks")
#'
#' # Compare with Shapiro-Wilk
#' env_sw <- TransformData(env_data,
#'                         data_type = "env",
#'                         test = "shapiro")
#'
#' # Example 10: Complete workflow
#'
#' # 1. Transform environmental variables
#' env_final <- PreData(TransformData(env_data,
#'                                    data_type = "env",
#'                                    scale = TRUE), NCols = 1)
#'
#' # Multivariate analysis
#' pca_result <- prcomp(env_final, scale = FALSE)
#' plot(pca_result)
#' biplot(pca_result)
#'
#' # 2. Transform biological data
#' bio_final <- PreData(TransformData(bio_data,
#'                                    data_type = "bio",
#'                                    trans = "r4"), NCols = 1)
#'
#' # Multivariate analysis
#' pca_result2 <- prcomp(bio_final, scale = FALSE)
#' plot(pca_result2)
#' biplot(pca_result2)
#'
#' @importFrom stats shapiro.test ks.test sd
#' @export
TransformData <- function(DF, shape = "w", data_type = "env", version = "l",
                          alpha = 0.05, min_sd = 0, test = "shapiro",
                          trans = "r4", scale = FALSE, rel_type = "total") {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'.\nIf your species are in columns, use shape = 'w'.\nIf your species are in the first column, use shape = 'l'")
  }

  if (!data_type %in% c("env", "bio", "percent", "pa", "scale", "rel")) {  # ← CORREGIDO
    stop("data_type must be 'env', 'bio', 'percent', 'pa', 'scale', or 'rel'")
  }

  if (!version %in% c("l", "s")) {
    stop("version must be 'l' or 's'")
  }

  if (!trans %in% c("r4", "l1", "sq", "Ch")) {
    stop("trans must be 'r4' (fourth root), 'l1' (log+1), 'sq' (sqrt), or 'Ch' (chord)")
  }

  if (!rel_type %in% c("total", "max")) {  # ← AGREGAR ESTA LÍNEA
    stop("rel_type must be 'total' or 'max'")
  }

  if (ncol(DF) < 2) {
    stop("DF must have at least 2 columns (identifiers + variables)")
  }

  # Data preparation based on format
  if (shape == "l") {
    # Long format: Species/variables in first column, stations in other columns
    Species <- DF[, 1]
    DF_numeric <- as.data.frame(t(DF[, -1]))
    colnames(DF_numeric) <- Species
    rownames(DF_numeric) <- colnames(DF[, -1])
    Stations <- rownames(DF_numeric)

    # Store first column for version "l" output
    first_col_name <- "Station"
    first_col_values <- Stations

  } else {
    # Wide format: Stations in first column, species/variables in other columns
    Stations <- DF[, 1]
    DF_numeric <- DF[, -1, drop = FALSE]
    Species <- colnames(DF_numeric)
    rownames(DF_numeric) <- Stations

    # Store first column for version "l" output
    first_col_name <- colnames(DF)[1]
    first_col_values <- Stations
  }

  # Select only numeric columns
  numeric_cols <- sapply(DF_numeric, is.numeric)
  if (sum(numeric_cols) == 0) {
    stop("DF must contain at least one numeric column")
  }
  DF_numeric <- DF_numeric[, numeric_cols, drop = FALSE]

  # Store original data for version "l" output
  original_numeric <- DF_numeric

  # ===============================================================
  # SCALE ONLY (no transformation)
  # ===============================================================
  if (data_type == "scale") {
    scaled_data <- apply_scale(DF_numeric)

    transformation_log <- data.frame(
      Variable = colnames(scaled_data),
      Transformation = "scale (z-score)",
      stringsAsFactors = FALSE
    )

    summary_stats <- data.frame(
      Metric = c("Input variables", "Scaled variables"),
      Count = c(ncol(DF_numeric), ncol(scaled_data))
    )

    if (version == "l") {
      output <- data.frame(
        temp_col = first_col_values,
        scaled_data,
        stringsAsFactors = FALSE
      )
      colnames(output)[1] <- first_col_name
      rownames(output) <- NULL

      attr(output, "transformation_report") <- transformation_log
      attr(output, "summary") <- summary_stats
      attr(output, "scaled") <- TRUE

      cat("\n")
      cat("========================================\n")
      cat("       SCALING SUMMARY (Z-SCORE)        \n")
      cat("========================================\n")
      print(summary_stats)
      cat("\n")

      return(output)
    } else {
      return(list(
        scaled_data = scaled_data,
        transformation_report = transformation_log,
        summary = summary_stats
      ))
    }
  }

  # ===============================================================
  # PRESENCE-ABSENCE
  # ===============================================================
  if (data_type == "pa") {
    pa_data <- apply_pa(DF_numeric)

    transformation_log <- data.frame(
      Variable = colnames(pa_data),
      Transformation = "presence-absence",
      stringsAsFactors = FALSE
    )

    summary_stats <- data.frame(
      Metric = c("Input variables", "Converted to PA"),
      Count = c(ncol(DF_numeric), ncol(pa_data))
    )

    if (version == "l") {
      output <- data.frame(
        temp_col = first_col_values,
        pa_data,
        stringsAsFactors = FALSE
      )
      colnames(output)[1] <- first_col_name
      rownames(output) <- NULL

      attr(output, "transformation_report") <- transformation_log
      attr(output, "summary") <- summary_stats

      cat("\n")
      cat("========================================\n")
      cat("   TRANSFORMATION SUMMARY (PA)           \n")
      cat("========================================\n")
      print(summary_stats)
      cat("\n")

      return(output)
    } else {
      return(list(
        pa_data = pa_data,
        transformation_report = transformation_log,
        summary = summary_stats
      ))
    }
  }

  # ===============================================================
  # PERCENTAGE/PROPORTION DATA
  # ===============================================================
  if (data_type == "percent") {

    # Convert percentages to proportions if needed
    DF_for_transform <- DF_numeric
    converted_vars <- character(0)

    for (var in colnames(DF_numeric)) {
      if (max(DF_numeric[[var]], na.rm = TRUE) > 1) {
        DF_for_transform[[var]] <- DF_numeric[[var]] / 100
        converted_vars <- c(converted_vars, var)
      }
    }

    if (length(converted_vars) > 0) {
      message("Note: Variables detected as percentages [0,100] converted to proportions [0,1]: ",
              paste(converted_vars, collapse = ", "))
    }

    # Validate all values are in range 0 and 1
    if (any(DF_for_transform < 0 | DF_for_transform > 1, na.rm = TRUE)) {
      stop("All values must be proportions [0,1] or percentages [0,100]")
    }

    # Apply arcsine-sqrt transformation
    arcsin_data <- apply_arcsin(DF_for_transform)

    # Apply scaling if requested
    if (scale) {
      arcsin_data <- apply_scale(arcsin_data)
    }

    trans_name <- if(scale) "arcsine-sqrt (scaled)" else "arcsine-sqrt"

    transformation_log <- data.frame(
      Variable = colnames(arcsin_data),
      Transformation = trans_name,
      stringsAsFactors = FALSE
    )

    summary_stats <- data.frame(
      Metric = c("Input variables", "Transformed variables"),
      Count = c(ncol(DF_numeric), ncol(arcsin_data))
    )

    if (version == "l") {
      output <- data.frame(
        temp_col = first_col_values,
        arcsin_data,
        stringsAsFactors = FALSE
      )
      colnames(output)[1] <- first_col_name
      rownames(output) <- NULL

      attr(output, "transformation_report") <- transformation_log
      attr(output, "summary") <- summary_stats
      attr(output, "scaled") <- scale

      cat("\n")
      cat("========================================\n")
      cat("   TRANSFORMATION SUMMARY (PERCENT)     \n")
      cat("========================================\n")
      print(summary_stats)
      cat("\nTransformation applied:", trans_name, "\n")
      if (length(converted_vars) > 0) {
        cat("Variables converted: ", paste(converted_vars, collapse = ", "), "\n")
      }
      cat("\n")

      return(output)
    } else {
      return(list(
        transformed_data = arcsin_data,
        transformation_report = transformation_log,
        summary = summary_stats
      ))
    }
  }

  # ===============================================================
  # RELATIVIZATION
  # ===============================================================
  if (data_type == "rel") {

    if (!rel_type %in% c("total", "max")) {
      stop("rel_type must be 'total' or 'max'")
    }

    # Apply relativization by row
    rel_data <- as.data.frame(t(apply(DF_numeric, 1, function(x) {
      if (rel_type == "total") {
        if (sum(x, na.rm = TRUE) == 0) return(x)
        x / sum(x, na.rm = TRUE)
      } else {
        if (max(x, na.rm = TRUE) == 0) return(x)
        x / max(x, na.rm = TRUE)
      }
    })))

    trans_name <- paste0("relativize_", rel_type)

    transformation_log <- data.frame(
      Variable = colnames(DF_numeric),
      Transformation = trans_name,
      stringsAsFactors = FALSE
    )

    summary_stats <- data.frame(
      Metric = c("Input variables", "Relativized variables"),
      Count = c(ncol(DF_numeric), ncol(rel_data))
    )

    if (version == "l") {
      output <- data.frame(
        temp_col = first_col_values,
        rel_data,
        stringsAsFactors = FALSE
      )
      colnames(output)[1] <- first_col_name
      rownames(output) <- NULL

      attr(output, "transformation_report") <- transformation_log
      attr(output, "summary") <- summary_stats
      attr(output, "rel_type") <- rel_type

      cat("\n")
      cat("========================================\n")
      cat("     RELATIVIZATION SUMMARY             \n")
      cat("========================================\n")
      print(summary_stats)
      cat("\nRelativization type:", rel_type, "\n\n")

      return(output)
    } else {
      return(list(
        transformed_data = rel_data,
        transformation_report = transformation_log,
        summary = summary_stats
      ))
    }
  }

  # ===============================================================
  # BIOLOGICAL DATA
  # ===============================================================
  if (data_type == "bio") {

    trans_func <- switch(trans,
                         "r4" = apply_fourth,
                         "l1" = apply_log,
                         "sq" = apply_sqrt,
                         "Ch" = apply_chord)

    trans_name <- switch(trans,
                         "r4" = "fourth_root",
                         "l1" = "log_p1",
                         "sq" = "sqrt",
                         "Ch" = "chord")

    transformed_data <- trans_func(DF_numeric)

    if (scale) {
      transformed_data <- apply_scale(transformed_data)
      trans_name <- paste0(trans_name, " (scaled)")
    }

    transformation_log <- data.frame(
      Variable = colnames(DF_numeric),
      Transformation = trans_name,
      stringsAsFactors = FALSE
    )

    summary_stats <- data.frame(
      Metric = c("Input variables", "Transformed variables"),
      Count = c(ncol(DF_numeric), ncol(transformed_data))
    )

    if (version == "l") {
      output <- data.frame(
        temp_col = first_col_values,
        transformed_data,
        stringsAsFactors = FALSE
      )
      colnames(output)[1] <- first_col_name
      rownames(output) <- NULL

      attr(output, "transformation_report") <- transformation_log
      attr(output, "summary") <- summary_stats
      attr(output, "scaled") <- scale

      cat("\n")
      cat("========================================\n")
      cat("   TRANSFORMATION SUMMARY (BIO)         \n")
      cat("========================================\n")
      print(summary_stats)
      cat("\nTransformation applied:", trans_name, "\n\n")

      return(output)
    } else {
      return(list(
        transformed_data = transformed_data,
        transformation_report = transformation_log,
        summary = summary_stats
      ))
    }
  }

  # ===============================================================
  # ENVIRONMENTAL DATA
  # ===============================================================
  if (data_type == "env") {

    # Step 1: Filter by standard deviation
    sds <- sapply(DF_numeric, sd, na.rm = TRUE)
    keep_vars <- names(sds)[sds > min_sd]

    if (length(keep_vars) == 0) {
      stop("No variables with SD > min_sd")
    }

    DF_filtered <- DF_numeric[, keep_vars, drop = FALSE]
    excluded_low_sd <- setdiff(names(DF_numeric), keep_vars)

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

    normal_dfs <- list()
    current_non_normal <- DF_filtered

    # Step 2: Initial normality test
    initial_results <- test_all_variables(current_non_normal, test)

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

    if (ncol(current_non_normal) > 0) {
      transformation_log$Initial_p[transformation_log$Variable %in% colnames(current_non_normal)] <-
        initial_results$p_values[colnames(current_non_normal)]
    }

    # Transformation sequence
    transform_sequence <- list(
      list(name = "sqrt", func = apply_sqrt),
      list(name = "log", func = apply_log),
      list(name = "fourth", func = apply_fourth)
    )

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

    non_normal_vars <- character(0)
    if (ncol(current_non_normal) > 0) {
      remaining_vars <- colnames(current_non_normal)
      transformation_log$Status[transformation_log$Variable %in% remaining_vars] <- "Non-normal"
      final_test <- test_all_variables(current_non_normal, test)
      transformation_log$Final_p[transformation_log$Variable %in% remaining_vars] <-
        final_test$p_values[remaining_vars]
      non_normal_vars <- remaining_vars
    }

    # Combine normal data
    if (length(normal_dfs) > 0) {
      normal_data <- do.call(cbind, normal_dfs)
    } else {
      normal_data <- data.frame(matrix(ncol = 0, nrow = nrow(DF_filtered)))
      rownames(normal_data) <- rownames(DF_filtered)
    }

    # Non-normal data (original values)
    if (length(non_normal_vars) > 0) {
      non_normal_data <- original_data[, non_normal_vars, drop = FALSE]
    } else {
      non_normal_data <- data.frame(matrix(ncol = 0, nrow = nrow(DF_filtered)))
      rownames(non_normal_data) <- rownames(DF_filtered)
    }

    # Combine all data
    if (ncol(normal_data) > 0 && ncol(non_normal_data) > 0) {
      combined_data <- cbind(normal_data, non_normal_data)
    } else if (ncol(normal_data) > 0) {
      combined_data <- normal_data
    } else if (ncol(non_normal_data) > 0) {
      combined_data <- non_normal_data
    } else {
      combined_data <- data.frame(row.names = rownames(DF_filtered))
    }

    # Apply scaling if requested
    if (scale) {
      combined_data <- apply_scale(combined_data)
    }

    transformation_log$Initial_p <- round(transformation_log$Initial_p, 4)
    transformation_log$Final_p <- round(transformation_log$Final_p, 4)
    transformation_log <- transformation_log[order(transformation_log$Status, decreasing = TRUE), ]
    rownames(transformation_log) <- NULL

    transform_counts <- table(transformation_log$Transformation[transformation_log$Status == "Normal"])

    summary_stats <- data.frame(
      Metric = c("Input variables",
                 "Excluded (low SD)",
                 "Initially normal",
                 names(transform_counts),
                 "Remaining non-normal",
                 "Total normalized",
                 if(scale) "Scaled (all variables)" else NULL),
      Count = c(
        length(all_vars) + length(excluded_low_sd),
        length(excluded_low_sd),
        sum(transformation_log$Transformation == "none" & transformation_log$Status == "Normal"),
        as.numeric(transform_counts),
        length(non_normal_vars),
        ncol(normal_data),
        if(scale) ncol(combined_data) else NULL
      )
    )

    if (version == "l") {
      output <- data.frame(
        temp_col = first_col_values,
        combined_data,
        stringsAsFactors = FALSE
      )
      colnames(output)[1] <- first_col_name
      rownames(output) <- NULL

      attr(output, "transformation_report") <- transformation_log
      attr(output, "summary") <- summary_stats
      attr(output, "normal_vars") <- colnames(normal_data)
      attr(output, "non_normal_vars") <- non_normal_vars
      attr(output, "scaled") <- scale

      cat("\n")
      cat("========================================\n")
      cat("     TRANSFORMATION SUMMARY (ENV)       \n")
      cat("========================================\n")
      print(summary_stats)
      cat("\nTransformation Report:\n")
      print(transformation_log)
      if (scale) {
        cat("\nNote: All variables (normal and non-normal) were scaled (z-score).\n")
      }
      cat("\n")

      return(output)
    } else {
      return(list(
        normal_data = if(scale) combined_data[, colnames(normal_data), drop = FALSE] else normal_data,
        non_normal_data = if(scale) combined_data[, non_normal_vars, drop = FALSE] else non_normal_data,
        transformation_report = transformation_log,
        summary = summary_stats,
        normal_vars = colnames(normal_data),
        non_normal_vars = non_normal_vars
      ))
    }
  }
}


# Helper functions (unchanged)
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


