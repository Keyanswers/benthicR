# benthicR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/yourusername/benthicR)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**benthicR** provides comprehensive tools for analyzing benthic marine and freshwater communities, including diversity indices, functional trait analysis, biomass conversions, and publication-ready visualizations.

## Installation

You can install the development version of benthicR from GitHub:
```r
# install.packages("devtools")
devtools::install_github("yourusername/benthicR")
```

## Features

### 🔬 Diversity Indices
Calculate taxonomic and functional diversity:
- Shannon, Simpson, Richness, Evenness
- Beta diversity (Jaccard, Sørensen, Bray-Curtis)
- LCBD (Local Contribution to Beta Diversity)
- SCBD (Species Contribution to Beta Diversity)

### 🧬 Functional Traits
Analyze functional diversity using fuzzy-coded trait matrices:
- Aggregate species by functional traits
- Calculate functional group abundances
- Support for multiple trait categories (feeding, mobility, body size, etc.)

### ⚖️ Biomass Conversions
Convert between biomass units with taxon-specific factors:
- Wet mass ↔ Dry mass ↔ AFDM
- Biomass to energy content (Joules)
- Built-in database of 47 taxonomic groups
- Support for custom conversion factors

### 📊 Visualization
Create publication-ready plots:
- Radar charts (single and multi-variable)
- Community composition pie charts and bar plots
- Interactive plots with plotly
- Custom color palettes for major taxonomic groups

### 🗺️ Mapping
Visualize sampling locations:
- Quick shapefile plotting with north arrow and scale bar
- Customizable map elements
- Included example dataset (Iceland administrative boundaries)

### 🧹 Data Preparation
Prepare and clean ecological data:
- Reshape abundance matrices
- Normality tests and data transformations
- Filter rare species by occurrence or abundance
- Multiple filtering methods (AND/OR logic)

## Quick Start
```r
library(benthicR)

# Example 1: Calculate diversity indices
data <- data.frame(
  Station = paste0("St", 1:5),
  Sp1 = c(10, 15, 8, 12, 20),
  Sp2 = c(5, 8, 12, 6, 10),
  Sp3 = c(20, 18, 15, 22, 25)
)

print(data)

rowSums(data[,2:3])

diversity <- Div(data, n =13, shape = "w")
print(diversity)

# Example 2: Explore community composition
community <- data.frame(
  Class = c("Polychaeta", "Bivalvia", "Gastropoda","Ophiuroidea"),
  Species = c("Capitella", "Nucula", "Nassarius", "Amphiura"),
  St1 = c(100, 70, 80, 44),
  St2 = c(80, 73, 85, 66)
)

print(community)

result <- ExploreComm(community, shape = "w", cutoff = 99, palette = MCol)
result
# Returns: summary tables, pie chart, and bar plot

# Example 3: Convert biomass units
biomass <- data.frame(
  Group = c("Amphipoda", "Bivalvia", "Polych. Errantia"),
  Species = c("Ampelisca", "Nucula", "Glycera"),
  St1_WM = c(100, 200, 150)  # Wet mass in mg
)

print(biomass)

# View available conversion factors
view_conversion_factors(conversion = "WM_to_DM")

# Convert wet mass to dry mass
biomass_DM <- ConvertBiomass(
  DF = biomass,
  taxon_col = "Group",
  biomass_cols = 3,
  conversion = "WM_to_DM"
)

print(biomass_DM)

# Example 4: Functional trait analysis
traits <- data.frame(
  Phylum = c("Annelida", "Mollusca", "Arthropoda"),
  Species = c("Capitella", "Nucula", "Ampelisca"),
  Predator = c(0.1, 0.0, 0.8),
  Sediment = c(0.9, 0.8, 0.2),
  Infauna = c(1.0, 1.0, 0.8)
)

print(traits)

abundance <- data.frame(
  Phylum = c("Annelida", "Mollusca", "Arthropoda"),
  Species = c("Capitella", "Nucula", "Ampelisca"),
  St1 = c(100, 50, 30),
  St2 = c(80, 60, 40)
)

print(abundance)

trait_abund <- Functional(
  Tr = traits,
  Dat = abundance,
  abundance_cols = 3:4,
  threshold = 0.1
)

print(trait_abund)

# Example 5: Create radar charts
diversity_data <- data.frame(
  Shannon = c(2.5, 2.3, 2.7, 3.1, 3.2),
  Richness = c(15, 18, 14, 23, 25),
  Evenness = c(0.8, 0.75, 0.85, 0.88, 0.91),
  row.names = c("Site A", "Site B", "Site C", "Site D", "Site E")
)

print(diversity_data)

RadChart(diversity_data, layout = c(1, 3), tsize = 2)

# Example 6: Map sampling stations
library(sf)
library(ggspatial)

iceland <- get_example_shapefile()
ShapeMap(iceland,
         xlim = c(-27.9, -12.5),
         ylim = c(63, 67))
```

## Color Palettes

benthicR includes built-in color palettes for common benthic taxonomic groups:
```r
# Marine groups
MCol()  # Mollusca, Annelida, Arthropoda, Echinodermata, etc.

# Worms
WCol()  # Polychaeta families

# Feeding guilds
FFCol()  # Predators, deposit feeders, suspension feeders, etc.

# Mixed functional groups
MFCol()  # Combined taxonomy and feeding mode
```

## Modules

| Module | Functions | Description |
|--------|-----------|-------------|
| **Diversity** | 10 functions | Taxonomic and beta diversity indices |
| **Functional Traits** | 1 function | Fuzzy-coded trait analysis |
| **Biomass Conversion** | 3 functions | Unit conversions with taxon-specific factors |
| **Visualization** | 3 functions | Radar charts, community plots |
| **Data Preparation** | 4 functions | Reshaping, filtering, transformations |
| **Color Palettes** | 4 functions | Taxonomic group colors |
| **Maps** | 2 functions | Shapefile visualization |

## Documentation

Access function documentation:
```r
?Div                    # Diversity indices
?Functional             # Functional trait analysis
?ConvertBiomass         # Biomass conversions
?RadChart               # Radar charts
?ExploreComm            # Community composition
?ShapeMap               # Map creation
?FilterRare             # Filter rare species
```

## Author

**Juan Carlos Rubio-Polania**  
📧 keyanswers@gmail.com

## Acknowledgments

Special thanks to:
- **ELROI** who lifts me up and  **Milo José**, my faithful four-legged companion during long coding sessions 🐕
- The R community for excellent packages that made benthicR possible
- All researchers who contributed conversion factors and ecological data to the literature

## License

MIT License - see [LICENSE](LICENSE) file for details

## Citation

If you use benthicR in your research, please cite:
```r
citation("benthicR")
```
```
Rubio-Polania, J.C. (2026). benthicR: Benthic Ecology Analysis Tools. 
R package version 0.1.0.
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue on GitHub.

## References

- Bache S, Wickham H (2022). *magrittr: A Forward-Pipe Operator for R*. R package version 2.0.3. https://CRAN.R-project.org/package=magrittr
- Baselga A, Orme D, Villeger S, De Bortoli J, Leprieur F, Logez M, Martinez-Santalla S, Martin-Devasa R, Gomez-Rodriguez C, Crujeiras R (2023). *betapart: Partitioning Beta Diversity into Turnover and Nestedness Components*. R package version 1.6. https://CRAN.R-project.org/package=betapart
- Chang W (2023). *extrafont: Tools for Using Fonts*. R package version 0.19. https://CRAN.R-project.org/package=extrafont
- Dray S, et al. (2024). *adespatial: Multivariate Multiscale Spatial Analysis*. R package version 0.3-24. https://CRAN.R-project.org/package=adespatial
- Dunnington D (2023). *ggspatial: Spatial Data Framework for ggplot2*. R package version 1.1.9. https://CRAN.R-project.org/package=ggspatial
- Garnier S, Ross N, Rudis R, Camargo AP, Sciaini M, Scherer C (2024). *viridis: Colorblind-Friendly Color Maps for R*. R package version 0.6.5.
- Guénard G, Legendre P (2022). Hierarchical Clustering with Contiguity Constraint in R. *Journal of Statistical Software*, 103(7), 1-26. https://doi.org/10.18637/jss.v103.i07
- Kuhn M (2008). Building Predictive Models in R Using the caret Package. *Journal of Statistical Software*, 28(5), 1-26. https://doi.org/10.18637/jss.v028.i05
- Nakazawa M (2024). *fmsb: Functions for Medical Statistics Book with some Demographic Data*. R package version 0.7.6. https://CRAN.R-project.org/package=fmsb
- Pebesma E (2018). Simple Features for R: Standardized Support for Spatial Vector Data. *The R Journal*, 10(1), 439-446. https://doi.org/10.32614/RJ-2018-009
- Pebesma E, Bivand R (2023). *Spatial Data Science: With Applications in R*. Chapman and Hall/CRC. https://doi.org/10.1201/9780429459016
- Sievert C (2020). *Interactive Web-Based Data Visualization with R, plotly, and shiny*. Chapman and Hall/CRC.
- Wickham H (2011). testthat: Get Started with Testing. *The R Journal*, 3, 5-10.
- Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.
- Wickham H, Pedersen T, Seidel D (2025). *scales: Scale Functions for Visualization*. R package version 1.4.0. https://CRAN.R-project.org/package=scales


## See Also

- [vegan](https://cran.r-project.org/package=vegan) - Community ecology analyses
- [betapart](https://cran.r-project.org/package=betapart) - Beta diversity partitioning
- [FD](https://cran.r-project.org/package=FD) - Functional diversity
- [marmap](https://cran.r-project.org/package=marmap) - Bathymetric data and oceanographic analysis
