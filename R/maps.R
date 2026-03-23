#' Get Path to Example Shapefile
#'
#' Returns the path to the example Iceland administrative boundaries shapefile
#' included with benthicR.
#'
#' @return Character string with path to example shapefile
#'
#' @examples
#' \dontrun{
#' # Get example shapefile path
#' iceland_shp <- get_example_shapefile()
#'
#' # Use with ShapeMap
#' ShapeMap(iceland_shp)
#' }
#'
#' @export
get_example_shapefile <- function() {
  shp_path <- system.file("extdata", "iceland_admin", "ISL_adm1.shp",
                          package = "benthicR")

  if (shp_path == "") {
    stop("Example shapefile not found. Please reinstall benthicR.")
  }

  return(shp_path)
}


#' Plot Shapefile with Map Elements
#'
#' Creates a map from a shapefile with optional coordinate limits, north arrow,
#' and scale bar. This is a wrapper around sf and ggplot2 for quick shapefile
#' visualization.
#'
#' @param shapefile_path Character. Full path to the shapefile (.shp file).
#' @param xlim Numeric vector of length 2. Longitude limits c(min, max).
#'   If NULL, uses full extent (default = NULL).
#' @param ylim Numeric vector of length 2. Latitude limits c(min, max).
#'   If NULL, uses full extent (default = NULL).
#' @param add_north Logical. Add north arrow? (default = TRUE).
#' @param add_scale Logical. Add scale bar? (default = TRUE).
#' @param north_pos Character. North arrow position: "tl" (top-left), "tr" (top-right),
#'   "bl" (bottom-left), "br" (bottom-right) (default = "tl").
#' @param scale_pos Character. Scale bar position: same options as north_pos (default = "bl").
#' @param fill_color Character. Fill color for polygons (default = "grey90").
#' @param border_color Character. Border color for polygons (default = "black").
#' @param xlab Character. X-axis label (default = "Longitude").
#' @param ylab Character. Y-axis label (default = "Latitude").
#'
#' @return A ggplot object.
#'
#' @details
#' This function reads a shapefile and creates a publication-ready map with
#' optional north arrow and scale bar. Requires the sf and ggspatial packages.
#'
#' The function does NOT install packages automatically. Make sure you have
#' installed sf and ggspatial before using:
#' \code{install.packages(c("sf", "ggspatial"))}
#'
#' @examples
#' \dontrun{
#' # Example 1: Use included Iceland shapefile
#' iceland_path <- get_example_shapefile()
#' ShapeMap(iceland_path)
#'
#' # Example 2: Zoom to specific region
#' ShapeMap(iceland_path,
#'          xlim = c(-27.9, -12.5),
#'          ylim = c(63, 67))
#'
#' # Example 3: Custom colors
#' ShapeMap(iceland_path,
#'          xlim = c(-27.9, -12.5),
#'          ylim = c(63, 67),
#'          fill_color = "lightblue",
#'          border_color = "navy",
#'          north_pos = "tr")
#'
#' # Example 4: Your own shapefile
#' ShapeMap("path/to/your/shapefile.shp")
#' }
#'
#' @seealso \code{\link[sf]{st_read}} for reading shapefiles,
#'   \code{\link[ggplot2]{geom_sf}} for plotting spatial data
#'
#' @importFrom ggplot2 ggplot geom_sf coord_sf labs theme_minimal unit
#' @export
ShapeMap <- function(shapefile_path,
                     xlim = NULL,
                     ylim = NULL,
                     add_north = TRUE,
                     add_scale = TRUE,
                     north_pos = "tl",
                     scale_pos = "bl",
                     fill_color = "grey90",
                     border_color = "black",
                     xlab = "Longitude",
                     ylab = "Latitude") {

  # Check if sf is installed
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required but not installed. Please install it with: install.packages('sf')")
  }

  # Check if ggspatial is installed (only if north arrow or scale needed)
  if ((add_north || add_scale) && !requireNamespace("ggspatial", quietly = TRUE)) {
    stop("Package 'ggspatial' is required for north arrow and scale bar. Please install it with: install.packages('ggspatial')")
  }

  # Validate shapefile path
  if (!file.exists(shapefile_path)) {
    stop("Shapefile not found at: ", shapefile_path)
  }

  # Read shapefile
  message("Reading shapefile: ", shapefile_path)
  shp <- sf::st_read(shapefile_path, quiet = TRUE)

  # Create base plot
  p <- ggplot2::ggplot(data = shp) +
    ggplot2::geom_sf(fill = fill_color, color = border_color) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal()

  # Add coordinate limits if provided
  if (!is.null(xlim) && !is.null(ylim)) {
    if (length(xlim) != 2 || length(ylim) != 2) {
      stop("xlim and ylim must be numeric vectors of length 2")
    }
    p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim)
  }

  # Add north arrow
  if (add_north) {
    p <- p + ggspatial::annotation_north_arrow(
      location = north_pos,
      which_north = "true",
      pad_x = ggplot2::unit(0.2, "cm"),
      pad_y = ggplot2::unit(0.2, "cm"),
      style = ggspatial::north_arrow_nautical,
      width = ggplot2::unit(2, "cm"),
      height = ggplot2::unit(2, "cm")
    )
  }

  # Add scale bar
  if (add_scale) {
    p <- p + ggspatial::annotation_scale(location = scale_pos)
  }

  return(p)
}
