#' @importFrom fmsb radarchart
#' @importFrom scales alpha
#' @importFrom caret preProcess
#' @importFrom plotly plot_ly layout
#' @importFrom ggplot2 ggplot aes geom_bar theme element_text coord_cartesian scale_y_continuous xlab ylab scale_fill_manual labs
#' @importFrom stats aggregate predict
#' @importFrom grDevices rainbow
#' @importFrom graphics legend par
#' @importFrom magrittr %>%
NULL

# Avoid R CMD check notes for ggplot2 aes() variables
utils::globalVariables(c("Species", "Relative_value", "Class"))

#' Exploratory Analysis of Community Composition
#'
#' Creates summary tables and visualizations (pie chart and bar plot) for
#' community composition data, showing taxonomic group abundances and
#' cumulative percentages. Colors are coordinated between pie chart and bar plot.
#'
#' @param DF A data.frame containing community data with different formats
#'   depending on the shape parameter:
#'   \itemize{
#'     \item If shape = "w" (wide): First column = taxonomic group/class,
#'       second column = species names, remaining columns = abundance data
#'       by station (stations in columns)
#'     \item If shape = "l" (long): First column = station names, second
#'       column = taxonomic group/class, third column = species names,
#'       remaining columns = abundance values (stations in rows)
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: taxonomic groups and species in rows,
#'       stations in columns (default)
#'     \item \code{"l"} - Long format: stations in rows, taxonomic groups
#'       and species identifiers in first columns
#'   }
#' @param label Character. Label for legend in bar plot (default = "Taxonomic Group").
#' @param cutoff Numeric. Cumulative percentage threshold for filtering species
#'   in bar plot (default = 95). Only species contributing to the first X% of
#'   cumulative abundance are shown.
#' @param palette Function. Color palette function from benthicR (MCol, WCol,
#'   FFCol, or MFCol). Default is MCol. These functions return named vectors
#'   mapping taxonomic groups to colors. If a group is not found in the palette,
#'   a default color will be assigned.
#'
#' @return A list with four elements:
#'   \itemize{
#'     \item \code{group_summary} - Total abundance by taxonomic group
#'     \item \code{species_summary} - Detailed table with species, abundances,
#'       relative abundance (\%), and cumulative abundance (\%)
#'     \item \code{pie_chart} - Interactive plotly pie chart of group composition
#'     \item \code{bar_plot} - ggplot2 bar plot of species composition
#'   }
#'
#' @details
#' This function performs exploratory analysis of community composition data,
#' providing both summary statistics and coordinated visualizations.
#'
#' \strong{Data format:}
#'
#' Wide format (shape = "w", default):
#' \preformatted{
#'   Class      Species           St1  St2  St3
#'   Annelida   Capitella         120   95  110
#'   Annelida   Glycera            45   38   52
#'   Mollusca   Nucula             23   18   31
#' }
#'
#' Long format (shape = "l"):
#' \preformatted{
#'   Station  Class      Species    Abundance
#'   St1      Annelida   Capitella  120
#'   St1      Annelida   Glycera     45
#'   St1      Mollusca   Nucula      23
#'   St2      Annelida   Capitella   95
#' }
#'
#' \strong{Color coordination:}
#'
#' Colors are obtained from the specified palette function (a named vector
#' dictionary mapping taxonomic groups to colors). These colors are then
#' consistently applied to both the pie chart and bar plot, ensuring that
#' each taxonomic group has the same color across all visualizations.
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item Processes data according to shape parameter
#'   \item Calculates total abundance per species across all stations
#'   \item Computes relative abundance (percentage of total)
#'   \item Calculates cumulative abundance (for dominance analysis)
#'   \item Maps colors from palette to taxonomic groups
#'   \item Creates pie chart showing taxonomic group composition
#'   \item Creates bar plot showing species contributing to cutoff\% of total abundance
#' }
#'
#' \strong{Cumulative abundance filtering:}
#'
#' The cutoff parameter filters species for the bar plot based on cumulative
#' abundance. For example, cutoff = 95 shows only species that together account
#' for 95\% of total abundance, highlighting the most important contributors
#' while reducing visual clutter from rare species.
#'
#' \strong{Available color palettes:}
#' \itemize{
#'   \item \code{MCol} - Marine color palette (default)
#'   \item \code{WCol} - Warm color palette
#'   \item \code{FFCol} - Full female color palette
#'   \item \code{MFCol} - Male-female color palette
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Wide format (default)
#' community_wide <- data.frame(
#'   Class = c("Annelida", "Annelida", "Mollusca", "Mollusca",
#'             "Arthropoda", "Arthropoda"),
#'   Species = c("Capitella capitata", "Glycera dibranchiata",
#'               "Nucula nucleus", "Macoma balthica",
#'               "Ampelisca brevicornis", "Corophium volutator"),
#'   St1 = c(120, 45, 23, 67, 34, 89),
#'   St2 = c(95, 38, 18, 54, 28, 72),
#'   St3 = c(110, 52, 31, 71, 41, 94)
#' )
#'
#' result <- ExploreComm(community_wide,
#'                       shape = "w",
#'                       label = "Taxonomic Class",
#'                       cutoff = 90)
#'
#' # View summary tables
#' result$group_summary
#' result$species_summary
#'
#' # View plots (colors match between pie and bar)
#' result$pie_chart    # Interactive
#' result$bar_plot     # Static
#'
#' # Example 2: Long format
#' community_long <- data.frame(
#'   Station = rep(c("St1", "St2", "St3"), each = 3),
#'   Class = rep(c("Annelida", "Mollusca", "Arthropoda"), 3),
#'   Species = rep(c("Capitella", "Nucula", "Ampelisca"), 3),
#'   Abundance = c(120, 23, 34, 95, 18, 28, 110, 31, 41)
#' )
#'
#' result_long <- ExploreComm(community_long,
#'                            shape = "l",
#'                            cutoff = 95)
#'
#' # Example 3: Using different color palettes
#' result_warm <- ExploreComm(community_wide,
#'                            shape = "w",
#'                            palette = WCol)
#'
#' result_female <- ExploreComm(community_wide,
#'                              shape = "w",
#'                              palette = FFCol)
#'
#' # Example 4: Show only top 80% of species
#' result_filtered <- ExploreComm(community_wide,
#'                                shape = "w",
#'                                cutoff = 80,
#'                                palette = MFCol)
#'
#' # Example 5: Save plots
#' result <- ExploreComm(community_wide, shape = "w", palette = MCol)
#'
#' # Save interactive pie chart
#' htmlwidgets::saveWidget(result$pie_chart, "composition_pie.html")
#'
#' # Save bar plot
#' ggsave("composition_bar.png", result$bar_plot,
#'        width = 12, height = 8, dpi = 300)
#' }
#'
#' @references
#' Clarke, K.R. & Warwick, R.M. (2001). Change in Marine Communities: An Approach
#' to Statistical Analysis and Interpretation (2nd ed.). PRIMER-E, Plymouth.
#'
#' @seealso \code{\link{Div}} for diversity indices,
#'   \code{\link{FilterRare}} for filtering rare species,
#'   \code{\link{MCol}}, \code{\link{WCol}}, \code{\link{FFCol}},
#'   \code{\link{MFCol}} for color palettes
#'
#' @importFrom plotly plot_ly layout
#' @importFrom ggplot2 ggplot aes geom_bar theme element_text coord_cartesian scale_y_continuous xlab ylab scale_fill_manual labs
#' @export
ExploreComm <- function(DF, shape = "w", label = "Taxonomic Group",
                        cutoff = 95, palette = MCol) {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'")
  }

  if (cutoff <= 0 || cutoff > 100) {
    stop("cutoff must be between 0 and 100")
  }

  if (!is.function(palette)) {
    stop("palette must be a color palette function (MCol, WCol, FFCol, or MFCol)")
  }

  # ============================================================================
  # DATA PREPARATION based on shape
  # ============================================================================

  if (shape == "w") {
    # Wide format: Class, Species, St1, St2, St3, ...
    if (ncol(DF) < 3) {
      stop("Wide format requires at least 3 columns: taxonomic group, species, and abundance data")
    }

    taxonomic_group <- DF[, 1]
    species_names <- DF[, 2]
    abundance_data <- DF[, -(1:2), drop = FALSE]

  } else {
    # Long format: Station, Class, Species, Abundance, ...
    if (ncol(DF) < 4) {
      stop("Long format requires at least 4 columns: station, taxonomic group, species, and abundance")
    }

    station_col <- DF[, 1]
    taxonomic_group <- DF[, 2]
    species_names <- DF[, 3]
    abundance_data <- DF[, -(1:3), drop = FALSE]

    # For long format, we need to aggregate by species across stations
    # Create a unique identifier for each species
    species_id <- paste(taxonomic_group, species_names, sep = "_")

    # Aggregate abundance by species
    agg_data <- aggregate(abundance_data,
                          by = list(Species_ID = species_id,
                                    Class = taxonomic_group,
                                    Species = species_names),
                          FUN = sum, na.rm = TRUE)

    # Reorder columns
    taxonomic_group <- agg_data$Class
    species_names <- agg_data$Species
    abundance_data <- agg_data[, -(1:3), drop = FALSE]
  }

  # ============================================================================
  # CALCULATE ABUNDANCES
  # ============================================================================

  # Calculate total abundance per species
  total_abundance <- rowSums(abundance_data, na.rm = TRUE)

  # Calculate grand total
  grand_total <- sum(total_abundance, na.rm = TRUE)

  # Create species summary table
  species_summary <- data.frame(
    Class = as.character(taxonomic_group),
    Species = as.character(species_names),
    Total = total_abundance,
    Relative_value = round((total_abundance / grand_total) * 100, 3),
    stringsAsFactors = FALSE
  )

  # Sort by relative abundance (descending)
  species_summary <- species_summary[order(-species_summary$Relative_value), ]

  # Calculate cumulative percentage
  species_summary$Cumulative_value <- cumsum(species_summary$Relative_value)

  # Reset row names
  rownames(species_summary) <- NULL

  # Group summary (for pie chart)
  group_summary <- aggregate(Relative_value ~ Class,
                             data = species_summary,
                             FUN = sum)
  colnames(group_summary) <- c("group", "value")

  # ============================================================================
  # COLOR MAPPING - Using palette dictionary
  # ============================================================================

  # Get the color dictionary from the palette function
  color_dict <- palette()

  # Get unique groups in the data
  unique_groups <- group_summary$group

  # Map colors to groups
  color_map <- character(length(unique_groups))
  names(color_map) <- unique_groups

  for (group in unique_groups) {
    if (group %in% names(color_dict)) {
      # Group found in palette
      color_map[group] <- color_dict[group]
    } else {
      # Group not in palette, assign a default color
      # Use a hash-based color generation for consistency
      color_map[group] <- rainbow(length(unique_groups))[which(unique_groups == group)]
    }
  }

  # ============================================================================
  # PIE CHART - Using mapped colors
  # ============================================================================

  # Colors for pie chart (in order of group_summary)
  pie_colors <- color_map[group_summary$group]

  pie_chart <- plotly::plot_ly(
    group_summary,
    labels = ~group,
    values = ~value,
    type = "pie",
    textposition = "inside",
    textinfo = "label+percent",
    insidetextfont = list(color = "#FFFFFF"),
    hoverinfo = "text",
    text = ~paste(group, "<br>", round(value, 2), "%"),
    textfont = list(family = "Arial", size = 16, color = "white"),
    marker = list(
      colors = pie_colors,
      line = list(color = "#FFFFFF", width = 1)
    )
  ) %>%
    plotly::layout(
      showlegend = TRUE,
      legend = list(font = list(family = "Arial", size = 14))
    )

  # ============================================================================
  # BAR PLOT - Using same colors as pie chart
  # ============================================================================

  # Filter species for bar plot based on cumulative threshold
  species_filtered <- species_summary[species_summary$Cumulative_value <= cutoff, ]

  # Sort by relative value for bar plot
  species_filtered <- species_filtered[order(species_filtered$Relative_value), ]

  # Make species a factor with ordered levels for proper plotting
  species_filtered$Species <- factor(
    species_filtered$Species,
    levels = species_filtered$Species
  )

  # Calculate y-axis max
  ymax <- max(species_filtered$Relative_value, na.rm = TRUE) + 2

  # Create bar plot with SAME colors as pie chart
  bar_plot <- ggplot2::ggplot(
    data = species_filtered,
    ggplot2::aes(x = Species, y = Relative_value, fill = Class)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14, family = "sans"),
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        face = "italic", colour = "black", size = 12
      ),
      axis.text.y = ggplot2::element_text(colour = "black", size = 12),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      aspect.ratio = 0.6
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::xlab("Species") +
    ggplot2::ylab("Relative abundance (%)") +
    ggplot2::scale_fill_manual(values = color_map) +  # SAME colors as pie
    ggplot2::labs(fill = label)

  # Return results
  result <- list(
    group_summary = group_summary,
    species_summary = species_summary,
    pie_chart = pie_chart,
    bar_plot = bar_plot
  )

  return(result)
}


#' Create Radar Charts for Multivariate Data
#'
#' Creates customizable radar charts (spider charts) for visualizing multivariate
#' data, with options for single or multiple charts, single or multiple variables
#' per chart, and extensive styling control.
#'
#' @param DF A data.frame or numeric vector containing the data to plot.
#'   Format depends on shape parameter:
#'   \itemize{
#'     \item If vector: creates a single radar chart (shape is ignored)
#'     \item If data.frame with shape = "w": stations in rows, variables in columns
#'     \item If data.frame with shape = "l": variables in rows, stations in columns
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations/samples in rows, variables in columns (default)
#'     \item \code{"l"} - Long format: variables in rows, stations/samples in columns
#'   }
#' @param Ftext Character. Font family to use for text (default = "Times New Roman").
#'   Note: Font must be available on your system. Use extrafont::fonts() to see
#'   available fonts.
#' @param tsize Numeric. Size of variable labels (default = 1.5).
#' @param col1 Character. Color for grid lines (default = "black").
#' @param col2 Character or vector. Line color for single series or when plotting
#'   multiple charts (default = "firebrick" for single, "dodgerblue" for multiple).
#' @param col3 Character vector. Colors for multiple series when MV is used.
#'   If not provided, uses default R colors (2:n) with transparency.
#' @param Vec Character. Title for single variable radar chart.
#' @param MV Character. Title for multi-variable radar chart.
#' @param layout Numeric vector of length 2. Layout for multiple plots as c(rows, cols).
#'   For example, c(1, 5) creates 1 row with 5 columns. Only used when creating
#'   multiple separate charts.
#' @param normalize Logical. Should data be normalized to 0-1 range when plotting
#'   multiple series? (default = TRUE when MV is provided). Normalization uses
#'   min-max scaling to make variables with different scales comparable.
#'
#' @return Invisibly returns NULL. The function is called for its side effect
#'   of creating plots in the active graphics device.
#'
#' @details
#' This function creates radar charts with three main usage modes:
#'
#' \strong{Mode 1: Single variable, single chart}
#'
#' When DF is a vector or when Vec parameter is provided, creates a single
#' radar chart showing one variable across multiple categories.
#'
#' \strong{Mode 2: Multiple variables, separate charts}
#'
#' When DF is a data.frame and neither MV nor Vec are provided, creates
#' separate radar charts for each variable (column in wide format, row in long format),
#' comparing values across samples.
#'
#' \strong{Mode 3: Multiple variables/samples, single chart}
#'
#' When MV parameter is provided, creates a single radar chart with multiple
#' overlaid series. Data are automatically normalized to 0-1 range for
#' comparability. Each row (wide format) or column (long format) becomes a
#' separate series in different colors.
#'
#' \strong{Data format (shape parameter):}
#'
#' Wide format (shape = "w", default):
#' \preformatted{
#'              Shannon  Richness  Evenness
#'   Station1   2.5      15        0.8
#'   Station2   2.3      18        0.75
#'   Station3   2.7      14        0.85
#' }
#'
#' Long format (shape = "l"):
#' \preformatted{
#'              Station1  Station2  Station3
#'   Shannon    2.5       2.3       2.7
#'   Richness   15        18        14
#'   Evenness   0.8       0.75      0.85
#' }
#'
#' \strong{Normalization:}
#'
#' In multi-series mode (MV provided), data are normalized using min-max scaling:
#' \deqn{x_{norm} = \frac{x - x_{min}}{x_{max} - x_{min}}}
#'
#' This ensures all variables are on the same 0-1 scale, making visual
#' comparison meaningful even when variables have different units or ranges.
#'
#' @examples
#' \dontrun{
#' # Example 1: Single variable radar chart
#' diversity_indices <- c(Shannon = 2.5, Richness = 15, Evenness = 0.8,
#'                        LCBD = 0.05, Dominance = 0.3)
#' RadChart(diversity_indices, Vec = "Station A Diversity")
#'
#' # Example 2: Multiple separate charts - wide format
#' diversity_wide <- data.frame(
#'   Shannon = c(2.5, 2.3, 2.7, 2.4, 2.6),
#'   Richness = c(15, 18, 14, 16, 17),
#'   Evenness = c(0.8, 0.75, 0.85, 0.78, 0.82),
#'   row.names = paste0("St", 1:5)
#' )
#'
#' RadChart(diversity_wide, shape = "w", layout = c(1, 3), col2 = "forestgreen")
#'
#' # Example 3: Multiple separate charts - long format
#' diversity_long <- data.frame(
#'   St1 = c(2.5, 15, 0.8),
#'   St2 = c(2.3, 18, 0.75),
#'   St3 = c(2.7, 14, 0.85),
#'   row.names = c("Shannon1", "Shannon2", "Shannon3")
#' )
#'
#' RadChart(diversity_long, shape = "l", layout = c(1, 3))
#'
#' # Example 4: Multiple series comparison - wide format
#' temporal_comparison <- data.frame(
#'   Winter = c(2.1, 2.3, 2.0),
#'   Spring = c(2.5, 2.7, 2.4),
#'   Summer = c(2.8, 3.0, 2.7),
#'   row.names = c("Station_A", "Station_B", "Station_C")
#' )
#'
#' RadChart(temporal_comparison,
#'          shape = "w",
#'          MV = "Shannon Diversity - Seasonal Variation",
#'          col3 = c("dodgerblue", "forestgreen", "firebrick"))
#'
#' # Example 5: Multiple series comparison - long format
#' temporal_long <- data.frame(
#'   Station_A = c(2.1, 2.5, 2.8),
#'   Station_B = c(2.3, 2.7, 3.0),
#'   Station_C = c(2.0, 2.4, 2.7),
#'   Station_D = c(2.5, 2.1, 3.1),
#'
#'   row.names = c("Winter", "Spring", "Summer")
#' )
#'
#' RadChart(temporal_long,
#'          shape = "l",
#'          MV = "Shannon Diversity - Seasonal",
#'          col3 = c("dodgerblue", "forestgreen", "firebrick"))
#' }
#'
#' @seealso \code{\link{Div}} for calculating diversity indices to plot,
#'   \code{\link[fmsb]{radarchart}} for the underlying plotting function,
#'   \code{\link[extrafont]{font_import}} for importing system fonts
#'
#' @importFrom fmsb radarchart
#' @importFrom scales alpha
#' @importFrom caret preProcess
#' @export
RadChart <- function(DF, shape = "w", Ftext = "Times New Roman", tsize = 1.5,
                     col1 = "black", col2, col3, Vec, MV,
                     layout, normalize = TRUE) {

  # Input validation
  if (!is.data.frame(DF) && !is.numeric(DF)) {
    stop("DF must be a data.frame or numeric vector")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'")
  }

  # Check and load fonts if available
  if (requireNamespace("extrafont", quietly = TRUE)) {
    tryCatch({
      if (length(extrafont::fonts()) > 0) {
        extrafont::loadfonts(device = "win", quiet = TRUE)
      }
    }, error = function(e) {
      message("Could not load custom fonts. Using default font.")
      Ftext <- ""
    })
  } else {
    Ftext <- ""
  }

  # ============================================================================
  # MODE 1: Single variable (vector input)
  # ============================================================================

  if (is.null(dim(DF))) {

    Max <- rep(max(DF, na.rm = TRUE), length(DF))
    Min <- rep(min(DF, na.rm = TRUE), length(DF))
    df1 <- data.frame(rbind(Max, Min, DF))

    cls1 <- round(seq(round(min(DF, na.rm = TRUE), digits = 2),
                      round(max(DF, na.rm = TRUE), digits = 2),
                      length.out = 5), digits = 2)

    if (missing(col2)) col2 <- "firebrick"
      if (missing(Vec)) Vec <- ""

    fn <- par(family = Ftext)
    on.exit(par(fn))

    fmsb::radarchart(df1,
                     title = Vec,
                     axislabcol = "black",
                     cglcol = col1,
                     pcol = col2,
                     pfcol = scales::alpha(col2, 0.4),
                     cglty = 7, cglwd = 1,
                     vlabels = colnames(df1),
                     vlcex = tsize,
                     plwd = 1.3,
                     calcex = 1.5,
                     axistype = 1,
                     caxislabels = cls1,
                     cex.main = 2.3)

    return(invisible(NULL))
  }

  # ============================================================================
  # TRANSPOSE DATA if needed (long format)
  # ============================================================================

  if (shape == "l") {
    DF <- as.data.frame(t(DF))
  }

  # From here on, DF is always in wide format (stations in rows, variables in columns)

  # ============================================================================
  # MODE 2: Multiple separate charts (one per variable)
  # ============================================================================

  if (missing(MV) && missing(Vec)) {

    # Calculate max and min for each variable
    Max <- apply(DF, 2, max, na.rm = TRUE)
    Min <- apply(DF, 2, min, na.rm = TRUE)

    # Create list of data frames (one per variable)
    DF_list <- list()
    for (i in seq_along(colnames(DF))) {
      DF_list[[i]] <- data.frame(
        rbind(
          rep(Max[i], nrow(DF)),
          rep(Min[i], nrow(DF)),
          DF[, i]
        )
      )
      colnames(DF_list[[i]]) <- rownames(DF)
      rownames(DF_list[[i]]) <- c("Max", "Min", colnames(DF)[i])
    }

    # Calculate axis labels for each variable
    cls_list <- list()
    for (i in colnames(DF)) {
      cls_list[[i]] <- round(seq(round(min(DF[, i], na.rm = TRUE), digits = 2),
                                 round(max(DF[, i], na.rm = TRUE), digits = 2),
                                 length.out = 5), digits = 2)
    }

    # Set up layout if provided
    if (!missing(layout)) {
      if (length(layout) != 2) {
        stop("layout must be a vector of length 2: c(rows, cols)")
      }
      par(mfrow = layout)
    }

    # Set margins
    par(mar = c(1.5, 0, 2.5, 0))

    if (missing(col2)) col2 <- "dodgerblue"

    fn <- par(family = Ftext)
    on.exit(par(fn))

    # Create one chart per variable
    for (i in seq_along(DF_list)) {
      fmsb::radarchart(DF_list[[i]],
                       title = rownames(DF_list[[i]])[3],
                       axislabcol = "black",
                       cglcol = col1,
                       pcol = col2,
                       pfcol = scales::alpha(col2, 0.5),
                       cglty = 7, cglwd = 1,
                       vlcex = tsize,
                       plwd = 1.3,
                       calcex = 1.5,
                       axistype = 1,
                       caxislabels = cls_list[[i]],
                       cex.main = 2.3)
    }

    return(invisible(NULL))
  }

  # ============================================================================
  # MODE 3: Multiple series in single chart (comparison mode)
  # ============================================================================

  # Normalize data if requested
  if (normalize) {
    R <- caret::preProcess(DF, method = "range")
    DF <- predict(R, DF)
    DF <- as.data.frame(DF)
  }

  # Calculate max and min
  Max <- apply(DF, 2, max, na.rm = TRUE)
  Min <- apply(DF, 2, min, na.rm = TRUE)

  DD <- rbind(Max, Min, DF)
  rownames(DD)[1:2] <- c("Max", "Min")

  cls <- round(seq(min(DD, na.rm = TRUE), max(DD, na.rm = TRUE),
                   length.out = 5), digits = 2)

  # Handle colors
  n_series <- nrow(DF)

  if (missing(col3)) {
    col3 <- 2:(n_series + 1)
  }

  # Apply transparency
  col3_fill <- scales::alpha(col3, 0.5)
  col3_line <- scales::alpha(col3, 1)

  if (missing(MV)) MV <- ""

  fn <- par(family = Ftext)
  on.exit(par(fn))

  fmsb::radarchart(DD,
                   title = MV,
                   axislabcol = "black",
                   pcol = col3_line,
                   cglcol = col1,
                   pfcol = col3_fill,
                   cglty = 7, cglwd = 2,
                   vlcex = tsize,
                   plty = 1,
                   plwd = 3,
                   calcex = 1.5,
                   axistype = 1,
                   caxislabels = cls,
                   cex.main = 2.3)

  # Add legend
  legend("topright",
         legend = rownames(DD)[-(1:2)],
         bty = "n",
         col = col3_line,
         pch = 19,
         cex = 2,
         pt.cex = 2)

  return(invisible(NULL))
}


#' Create Multiple Radar Charts with Multiple Series Each
#'
#' Creates multiple radar charts (one per index/variable group) where each chart
#' shows multiple series (time points, treatments, conditions) in different colors.
#' Each index gets its own independent scale (0 to maximum value of that index).
#' This is specialized for comparing the same index measured under different
#' conditions across multiple stations.
#'
#' @param DF A data.frame containing the data to plot.
#'   \itemize{
#'     \item If shape = "w": Stations in rows, variables in columns. Variables
#'       will be transposed to create radars where stations form the axes.
#'     \item If shape = "l": Variables in rows, stations in columns. First column
#'       can be named "sta" with variable names, or rownames contain variable names.
#'   }
#' @param shape Data format specification (default = "w"):
#'   \itemize{
#'     \item \code{"w"} - Wide format: stations in rows, variables in columns (default)
#'     \item \code{"l"} - Long format: variables in rows, stations in columns
#'   }
#' @param layout Numeric vector of length 2. Layout for multiple plots as c(rows, cols).
#'   If NULL, automatically calculates optimal layout (default = NULL).
#' @param col Character or character vector. Colors for the series:
#'   \itemize{
#'     \item \code{"rainbow"} - Uses rainbow color palette (default)
#'     \item \code{"viridis"} - Uses viridis color palette (requires viridis package)
#'     \item Character vector - Custom colors (e.g., c("red", "blue", "green"))
#'   }
#' @param Ftext Character. Font family to use for text (default = "Times New Roman").
#' @param tsize Numeric. Size of text elements including variable labels, title,
#'   and legend (default = 1.1).
#' @param col1 Character. Color for grid lines (default = "grey80").
#' @param add_legend Logical. Should a legend be added to each radar? (default = TRUE).
#' @param verbose Logical. Should group detection information be printed? (default = FALSE).
#'
#' @return Invisibly returns a list with:
#'   \itemize{
#'     \item \code{groups} - Names of index groups detected
#'     \item \code{data_per_group} - List of data frames, one per group
#'     \item \code{num_radars} - Number of radars created
#'     \item \code{layout_used} - Layout dimensions used c(rows, cols)
#'   }
#'
#' @details
#' This function is specialized for creating multiple radar charts where:
#' \itemize{
#'   \item Each radar represents one ecological index (Shannon, Richness, FRic, etc.)
#'   \item Within each radar, multiple series are shown in different colors
#'     (e.g., different time points, treatments, or data sources)
#'   \item Each index has its own independent scale: minimum = 0, maximum = highest
#'     value observed for that index
#'   \item Stations/samples form the axes of each radar
#' }
#'
#' \strong{Variable naming and grouping:}
#'
#' Variables are automatically grouped by extracting the base name (everything
#' before the first number or underscore):
#' \itemize{
#'   \item Shannon_Winter, Shannon_Spring → group "Shannon"
#'   \item Shannon1, Shannon2, Shannon3 → group "Shannon"
#'   \item Div_2020, Div_2021 → group "Div"
#'   \item FRic1, FRic2 → group "FRic"
#' }
#'
#' The extraction stops at the first digit (0-9) or underscore (_) encountered.
#'
#' \strong{Scale calculation:}
#'
#' Each radar uses a scale from 0 to the maximum value observed for that specific
#' index across all series. This ensures:
#' \itemize{
#'   \item Shannon (range 2-3) has appropriate spacing with scale 0-3
#'   \item Richness (range 10-25) has appropriate spacing with scale 0-25
#'   \item Different indices are not forced onto the same scale
#'   \item All radars start at 0 for consistent interpretation
#' }
#'
#' \strong{Data format examples:}
#'
#' Wide format (shape = "w", default) - stations in rows:
#' \preformatted{
#'      Shannon_Winter  Shannon_Spring  Richness_Winter  Richness_Spring
#' St1  2.1             2.5             15               18
#' St2  2.3             2.7             17               20
#' St3  2.0             2.4             14               17
#' }
#'
#' Long format (shape = "l") - variables in rows with "sta" column:
#' \preformatted{
#'   sta              St1   St2   St3   St4   St5
#'   Shannon_Winter   2.1   2.3   2.0   2.2   2.4
#'   Shannon_Spring   2.5   2.7   2.4   2.6   2.8
#'   Richness_Winter  15    17    14    16    18
#' }
#'
#' Long format (shape = "l") - variables in rows without "sta" column:
#' \preformatted{
#'                  St1   St2   St3   St4   St5
#' Shannon_Winter   2.1   2.3   2.0   2.2   2.4
#' Shannon_Spring   2.5   2.7   2.4   2.6   2.8
#' Richness_Winter  15    17    14    16    18
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Wide format - temporal diversity
#' diversity_temporal <- data.frame(
#'   Shannon_Winter = c(2.1, 2.3, 2.0, 2.2, 2.4),
#'   Shannon_Spring = c(2.5, 2.7, 2.4, 2.6, 2.8),
#'   Shannon_Summer = c(2.8, 3.0, 2.7, 2.9, 3.1),
#'   Richness_Winter = c(15, 17, 14, 16, 18),
#'   Richness_Spring = c(18, 20, 17, 19, 21),
#'   Richness_Summer = c(20, 22, 19, 21, 23),
#'   row.names = paste0("St", 1:5)
#' )
#'
#' RadChartMulti(diversity_temporal,
#'               shape = "w",
#'               layout = c(1, 2),
#'               col = c("dodgerblue", "forestgreen", "firebrick"))
#'
#' # Example 2: Wide format with numbered variables
#' diversity_numbered <- data.frame(
#'   Div1 = runif(10, 0, 5),
#'   Div2 = runif(10, 0, 6),
#'   Div3 = runif(10, 0, 4),
#'   Rich1 = sample(8:20, 10),
#'   Rich2 = sample(2:30, 10),
#'   Rich3 = sample(1:15, 10),
#'   row.names = letters[1:10]
#' )
#'
#' RadChartMulti(diversity_numbered,
#'               shape = "w",
#'               layout = c(1, 2),
#'               col = c("coral", "seagreen", "steelblue"))
#'
#' # Example 3: Long format with "sta" column
#' df_long <- data.frame(
#'   sta = c("Div1", "Div2", "Div3", "Rich1", "Rich2", "Rich3"),
#'   a = c(2.5, 3.1, 2.8, 15, 18, 16),
#'   b = c(2.3, 2.8, 2.6, 18, 20, 19),
#'   c = c(2.7, 3.3, 3.0, 14, 17, 15),
#'   d = c(2.4, 2.9, 2.7, 16, 19, 17)
#' )
#'
#' RadChartMulti(df_long,
#'               shape = "l",
#'               layout = c(1, 2),
#'               col = c("coral", "seagreen", "steelblue"),
#'               verbose = TRUE)  # Shows detected groups
#'
#' # Example 4: Long format without "sta" column (using rownames)
#' df_long2 <- data.frame(
#'   St1 = c(2.1, 2.5, 2.8, 15, 18, 20),
#'   St2 = c(2.3, 2.7, 3.0, 17, 20, 22),
#'   St3 = c(2.0, 2.4, 2.7, 14, 17, 19),
#'   row.names = c("Shannon_Winter", "Shannon_Spring", "Shannon_Summer",
#'                 "Richness_Winter", "Richness_Spring", "Richness_Summer")
#' )
#'
#' RadChartMulti(df_long2,
#'               shape = "l",
#'               layout = c(1, 2),
#'               col = "rainbow")
#'
#' # Example 5: Automatic layout with viridis colors
#' functional <- data.frame(
#'   FRic_2020 = c(25, 22, 24, 23),
#'   FRic_2021 = c(27, 24, 26, 25),
#'   FRic_2022 = c(26, 23, 25, 24),
#'   FEve_2020 = c(0.45, 0.48, 0.46, 0.47),
#'   FEve_2021 = c(0.48, 0.51, 0.49, 0.50),
#'   FEve_2022 = c(0.47, 0.50, 0.48, 0.49),
#'   FDiv_2020 = c(0.65, 0.68, 0.66, 0.67),
#'   FDiv_2021 = c(0.68, 0.71, 0.69, 0.70),
#'   FDiv_2022 = c(0.67, 0.70, 0.68, 0.69),
#'   row.names = c("North", "Central", "South", "East")
#' )
#'
#' RadChartMulti(functional,
#'               shape = "w",
#'               col = "viridis")  # Automatic 2x2 layout
#'
#' # Example 6: Custom colors and text size
#' restoration <- data.frame(
#'   Shannon_Before = c(1.8, 1.6, 1.9, 1.7),
#'   Shannon_After = c(2.5, 2.3, 2.6, 2.4),
#'   Richness_Before = c(12, 10, 13, 11),
#'   Richness_After = c(18, 16, 19, 17),
#'   Evenness_Before = c(0.65, 0.62, 0.67, 0.64),
#'   Evenness_After = c(0.78, 0.75, 0.80, 0.77),
#'   row.names = c("SiteA", "SiteB", "SiteC", "SiteD")
#' )
#'
#' RadChartMulti(restoration,
#'               layout = c(1, 3),
#'               col = c("firebrick", "forestgreen"),
#'               tsize = 1.3,
#'               col1 = "black")
#'
#' # Example 7: Treatment gradient
#' treatment_gradient <- data.frame(
#'   Richness_Control = c(18, 20, 19, 21, 17, 22),
#'   Richness_Low = c(16, 18, 17, 19, 15, 20),
#'   Richness_Medium = c(14, 16, 15, 17, 13, 18),
#'   Richness_High = c(12, 14, 13, 15, 11, 16),
#'   row.names = paste0("Rep", 1:6)
#' )
#'
#' RadChartMulti(treatment_gradient,
#'               layout = c(1, 1),
#'               col = c("forestgreen", "gold", "orange", "red"),
#'               Ftext = "Arial",
#'               tsize = 1.5)
#'
#' # Example 8: Retrieve and use output information
#' result <- RadChartMulti(diversity_temporal,
#'                         layout = c(1, 2),
#'                         verbose = TRUE)
#'
#' # See what groups were detected
#' print(result$groups)
#' # See layout used
#' print(result$layout_used)
#' # Access data for specific group
#' print(result$data_per_group$Shannon)
#' }
#'
#' @references
#' Saary, P., et al. (2016). RadarChart: Radar chart graphics. R package.
#'
#' @seealso \code{\link{RadChart}} for standard radar charts with more modes,
#'   \code{\link{Div}} for calculating diversity indices
#'
#' @importFrom fmsb radarchart
#' @importFrom scales alpha
#' @importFrom grDevices rainbow
#' @export
RadChartMulti <- function(DF, shape = "w", layout = NULL, col = "rainbow",
                          Ftext = "Times New Roman", tsize = 1.1,
                          col1 = "grey80", add_legend = TRUE, verbose = FALSE) {

  # Input validation
  if (!is.data.frame(DF)) {
    stop("DF must be a data.frame")
  }

  if (!shape %in% c("l", "w")) {
    stop("shape must be 'l' or 'w'")
  }

  # Check and load fonts if available
  if (requireNamespace("extrafont", quietly = TRUE)) {
    tryCatch({
      if (length(extrafont::fonts()) > 0) {
        extrafont::loadfonts(device = "win", quiet = TRUE)
      }
    }, error = function(e) {
      message("Could not load custom fonts. Using default font.")
      Ftext <- ""
    })
  } else {
    Ftext <- ""
  }

  # ============================================================================
  # HELPER FUNCTION: Extract base group name
  # ============================================================================

  extract_group_name <- function(name) {
    # Find position of first number or underscore
    # Pattern: any digit (0-9) or underscore (_)
    pos <- regexpr("[0-9_]", name)

    if (pos == -1) {
      # No number or underscore found, use complete name
      return(name)
    } else {
      # Extract everything before first match
      base_name <- substr(name, 1, pos - 1)
      return(base_name)
    }
  }

  # ============================================================================
  # PROCESS DATA BASED ON SHAPE
  # ============================================================================

  if (shape == "l") {
    # Long format: variables in rows, stations in columns
    if ("sta" %in% colnames(DF)) {
      # Has "sta" column with variable names
      dft <- data.frame(t(DF[, -1]))
      colnames(dft) <- DF[, 1]
    } else {
      # No "sta" column, assume rownames are variable names
      dft <- DF
    }

    # Identify groups by row names
    row_names <- rownames(dft)
    base_groups <- unique(sapply(row_names, extract_group_name))

    # Create list with data frames per group
    data_per_group <- list()
    for (group in base_groups) {
      # Find rows that start with base name
      group_rows <- grep(paste0("^", group), row_names, value = TRUE)
      if (length(group_rows) > 0) {
        data_per_group[[group]] <- dft[group_rows, , drop = FALSE]
      }
    }

  } else {
    # Wide format: stations in rows, variables in columns
    # Transpose to have variables in rows
    dft <- data.frame(t(DF))

    # Identify groups by row names
    row_names <- rownames(dft)
    base_groups <- unique(sapply(row_names, extract_group_name))

    # Create list with data frames per group
    data_per_group <- list()
    for (group in base_groups) {
      # Find rows that start with base name
      group_rows <- grep(paste0("^", group), row_names, value = TRUE)
      if (length(group_rows) > 0) {
        data_per_group[[group]] <- dft[group_rows, , drop = FALSE]
      }
    }
  }

  # Show detected groups information if verbose
  if (verbose) {
    message("\n=== DETECTED GROUPS ===")
    for (group in names(data_per_group)) {
      message("Group: ", group)
      message("  Rows: ", paste(rownames(data_per_group[[group]]), collapse = ", "))
    }
    message("========================\n")
  }

  # ============================================================================
  # SET UP COLORS
  # ============================================================================

  color_palettes <- list()

  if (is.character(col) && length(col) == 1) {
    if (col == "rainbow") {
      for (group in names(data_per_group)) {
        n_rows <- nrow(data_per_group[[group]])
        color_palettes[[group]] <- grDevices::rainbow(n_rows)
      }
    } else if (col == "viridis") {
      if (!requireNamespace("viridis", quietly = TRUE)) {
        stop("To use 'viridis' you need to install the viridis package: install.packages('viridis')")
      }
      for (group in names(data_per_group)) {
        n_rows <- nrow(data_per_group[[group]])
        color_palettes[[group]] <- viridis::viridis(n_rows)
      }
    } else {
      warning("col is a string but not 'rainbow' or 'viridis'. Using rainbow by default.")
      for (group in names(data_per_group)) {
        n_rows <- nrow(data_per_group[[group]])
        color_palettes[[group]] <- grDevices::rainbow(n_rows)
      }
    }
  } else if (is.vector(col) && length(col) > 0) {
    # Use provided color vector
    for (group in names(data_per_group)) {
      n_rows <- nrow(data_per_group[[group]])
      if (length(col) >= n_rows) {
        color_palettes[[group]] <- col[1:n_rows]
      } else {
        # Repeat colors if not enough
        color_palettes[[group]] <- rep(col, length.out = n_rows)
      }
    }
  } else {
    # Default: use rainbow
    for (group in names(data_per_group)) {
      n_rows <- nrow(data_per_group[[group]])
      color_palettes[[group]] <- grDevices::rainbow(n_rows)
    }
  }

  # ============================================================================
  # SET UP LAYOUT
  # ============================================================================

  num_plots <- length(data_per_group)

  if (is.null(layout)) {
    # Automatic layout
    rows <- ceiling(sqrt(num_plots))
    cols <- ceiling(num_plots / rows)
    if (verbose) {
      message("Automatic layout: ", rows, " rows x ", cols, " columns")
    }
  } else if (is.vector(layout) && length(layout) == 2) {
    rows <- layout[1]
    cols <- layout[2]

    total_space <- rows * cols
    if (total_space < num_plots) {
      warning("Layout ", rows, " x ", cols, " has ", total_space,
              " spaces but needs ", num_plots, " radars. Adjusting...")
      rows <- ceiling(sqrt(num_plots))
      cols <- ceiling(num_plots / rows)
      if (verbose) {
        message("Layout adjusted to: ", rows, " rows x ", cols, " columns")
      }
    }

    if (verbose) {
      message("Custom layout: ", rows, " rows x ", cols, " columns")
    }
  } else {
    stop("layout must be NULL or a vector c(rows, cols)")
  }

  par(mfrow = c(rows, cols), mar = c(2, 2, 4, 2))

  fn <- par(family = Ftext)
  on.exit(par(fn))

  # ============================================================================
  # CREATE RADARS
  # ============================================================================

  for (i in seq_along(data_per_group)) {
    group <- names(data_per_group)[i]
    group_data <- data_per_group[[group]]
    group_colors <- color_palettes[[group]]

    # Calculate max value for this group
    max_value <- max(group_data, na.rm = TRUE)

    # Create max/min rows: min always 0, max is the real maximum
    max_min <- data.frame(row.names = c("max", "min"))

    for (station in colnames(group_data)) {
      max_min[[station]] <- c(max_value, 0)
    }

    # Combine data
    plot_data <- rbind(max_min, group_data)

    # Create axis labels based on real maximum
    num_labels <- 5
    labels <- round(seq(0, max_value, length.out = num_labels), 1)

    # Create radar chart
    fmsb::radarchart(plot_data,
                     pcol = group_colors,
                     pfcol = scales::alpha(group_colors, 0.4),
                     plwd = 3,
                     plty = 1,
                     axistype = 1,
                     cglcol = col1,
                     cglty = 1,
                     axislabcol = "grey20",
                     caxislabels = labels,
                     vlcex = tsize,
                     title = group,
                     cex.main = tsize)

    # Add legend
    if (add_legend) {
      legend("topright",
             legend = rownames(group_data),
             bty = "n",
             pch = 20,
             col = group_colors,
             pt.cex = tsize,
             cex = tsize,
             horiz = FALSE)
    }
  }

  # ============================================================================
  # RETURN INFORMATION
  # ============================================================================

  return(invisible(list(
    groups = names(data_per_group),
    data_per_group = data_per_group,
    num_radars = num_plots,
    layout_used = c(rows, cols)
  )))
}
