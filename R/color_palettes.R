#' Color Palette for Marine Benthic Taxa
#'
#' Returns a named vector with predefined colors for the main taxonomic
#' groups of marine benthic communities.
#'
#' @return A named vector where names are taxa and values are color codes
#'   in hexadecimal format or R color names.
#'
#' @details
#' This palette is designed to maintain visual consistency in graphs
#' representing taxonomic composition of marine benthic communities.
#' Colors were selected to maximize contrast between groups and facilitate
#' visual interpretation.
#'
#' @examples
#' # Get the complete palette
#' marine_colors <- MCol()
#' marine_colors
#'
#' # Use in a plot
#' library(ggplot2)
#' data <- data.frame(
#'   taxa = c("Annelida", "Mollusca", "Crustacea"),
#'   abundance = c(45, 30, 25)
#' )
#' ggplot(data, aes(x = taxa, y = abundance, fill = taxa)) +
#'   geom_bar(stat = "identity") +
#'   scale_fill_manual(values = MCol())
#'
#' @export

MCol<-function()
{
  c(
    "Porifera" = "mediumpurple4",
    "Calcarea" = "mediumpurple3",
    "Demospongiae" = "mediumpurple2",
    "Hexactinellida" = "mediumpurple1",
    "Homoscleromorpha" = "thistle2",

    "Cnidaria" = "sienna4",
    "Anthozoa" = "sienna3",
    "Hydrozoa" = "peru",
    "Scyphozoa" = "burlywood3",


    "Platyhelminthes" = "goldenrod4",
    "Turbellaria" = "goldenrod2",

    "Nemertea" = "tan4",
    "Anopla" = "tan3",
    "Enopla" = "tan1",

    "Phoronida" = "khaki4",
    "Phoronis" = "khaki3",
    "Phoronopsis" = "khaki1",

    "Bryozoa" = "tan4",
    "Gymnolaemata" = "tan3",
    "Stenolaemata" = "tan2",

    "Brachiopoda" = "tan4",
    "Rhynchonellata" = "tan3",
    "Linguliformea" = "tan2",
    "Craniiformea" = "tan1",

    "Annelida" = "olivedrab4",
    "Polychaeta" = "olivedrab3",

    "Mollusca" = "seagreen4",
    "Gastropoda" = "seagreen3",
    "Bivalvia" = "darkseagreen3",
    "Polyplacophora" = "darkseagreen2",
    "Scaphopoda" = "darkseagreen1",
    "Cephalopoda" = "honeydew3",

    "Arthropoda" = "royalblue4",
    "Crustacea" = "royalblue3",
    "Malacostraca" = "royalblue2",
    "Decapoda" = "royalblue1",
    "Amphipoda" = "lightsteelblue3",
    "Isopoda" = "lightsteelblue2",
    "Stomatopoda" = "lightsteelblue1",

    "Thecostraca" = "steelblue3",
    "Copepoda" = "steelblue2",
    "Ostracoda" = "steelblue1",
    "Branchiopoda" = "lightblue",

    "Echinodermata" = "darkorange3",
    "Asteroidea" = "darkorange2",
    "Ophiuroidea" = "orange",
    "Echinoidea" = "orange2",
    "Holothuroidea" = "orange1",
    "Crinoidea" = "moccasin",

    "Chordata" = "indianred4",
    "Chondrichthyes" = "indianred3",
    "Osteichthyes" = "indianred2",

    "Ascidiacea" = "mistyrose3",
    "Cephalochordata" = "mistyrose2",

    "Total" = "firebrick4"

      )
}

#' Color Palette for Continental Benthic Taxa
#'
#' Returns a named vector with predefined colors for taxonomic groups
#' of benthic communities from continental systems (freshwater).
#'
#' @return A named vector where names are taxa and values are color codes.
#'
#' @details
#' This palette includes common taxonomic groups in freshwater and estuarine
#' environments, with colors selected to maximize visual contrast.
#'
#' @examples
#' # Get the palette
#' continental_colors <- WCol()
#'
#' # View all available groups
#' names(WCol())
#'
#' @export

WCol <-function ()
{
  c(

    "Porifera"      = "mediumpurple4",
    "Demospongiae"  = "mediumpurple2",

    "Cnidaria"      = "sienna4",
    "Hydrozoa"      = "peru",

    "Platyhelminthes" = "goldenrod4",

    "Turbellaria"     = "goldenrod2",

    "Nemertea" = "tan4",
    "Enopla"   = "tan2",

    "Nematoda"       = "burlywood4",
    "Adenophorea"    = "burlywood3",
    "Secernentea"    = "burlywood1",

    "Rotifera"     = "khaki4",
    "Bdelloidea"   = "khaki3",
    "Monogononta"  = "khaki1",

    "Gastrotricha"  = "darkkhaki",
    "Chaetonotida"  = "khaki2",

    "Tardigrada"    = "darkolivegreen4",
    "Eutardigrada"  = "darkolivegreen2",

    "Annelida"    = "olivedrab4",
    "Clitellata" = "olivedrab3",
    "Oligochaeta" = "olivedrab2",
    "Hirudinea"   = "olivedrab",


    "Mollusca"   = "seagreen4",
    "Gastropoda" = "seagreen3",
    "Bivalvia"   = "seagreen1",
    "Arthropoda" = "steelblue4",
    "Crustacea"  = "steelblue3",
    "Insecta"    = "steelblue2",
    "Arachnida"  = "steelblue1",

    "Total" = "firebrick4"
  )

}

#' Color Palette for Functional Traits
#'
#' Returns a named vector with predefined colors for different functional
#' trait categories in benthic communities.
#'
#' @return A named vector with color codes for functional traits related to
#'   feeding habit, mobility, position, habitat, volume, size, longevity,
#'   flexibility and protection.
#'
#'@details
#' Traits are organized by categories:
#' \itemize{
#'   \item Feeding: Pr (Predator), Om (omnivore), He (herbivore),
#'     Par (suspensiver-filter feeder), PP (primary producer), Sca (Scavenger), De (Detritivore),
#'     Sedi (sedimentivore-deposit feeder)
#'    \item Mobility: HMot (High mobility), Mot (Mobile), Sed (sedentary),
#'    Ses (sessile), Bur (burrower)
#'   \item Position: Epi (epifauna), Bi (In both), In (infauna)
#'   \item Habitat: Sur (surpercifial), Attached (Atta), Bu (burrow)
#'   \item Volume: VoSm (small volume), VoMe (medium volume), VoLa (large volume)
#'   \item Size: VLa (very large), La (large), Me (medium), Sm (small)
#'   \item Longevity: ALe5 (<5yr), B51 (5<yr<10), B10to50 (10<yr<50), Old50 (>50yr)
#'   \item Flexibility: Hi45 (Fl>45°), Lo (10°<Fl<45°), No (Fl<10°)
#'   \item Protection: HaVe (hard/heavy), Stro (strong), Npro (no protection),
#'     Fra (fragile)
#'}
#'
#' @examples
#' # Get trait palette
#' trait_colors <- MFCol()
#'
#' # View available categories
#' names(MFCol())
#'
#' @export

MFCol <-function () {
  c(
    "Pr"   = "#483D8B",
    "Om"   = "#104E8B",
    "He"   = "#00688B",
    "Par"  = "#1874CD",
    "PP"   = "#1C86EE",
    "Sca"  = "#1E90FF",
    "De"   = "#00B2EE",
    "Sedi" = "#009ACD",

    "HMot" = "#8B1A1A",
    "Mot"  = "#B22222",
    "Sed"  = "#CD2626",
    "Ses"  = "#EE2C2C",
    "Bur" = "#FF3030",

    "Epi"  = "#8B4500",
    "Bi"   = "#8B7500",
    "In"   = "#FF7F00",

    "Sur"  = "#8B6969",
    "Atta" = "#CD9B9B",
    "Bu"  = "#FFC1C1",

    "VoSm" = "#8B7500",
    "VoMe" = "#CDAD00",
    "VoLa" = "#EEC900",

    "VLa"  = "#5D478B",
    "La"   = "#7B68EE",
    "Med"  = "#0000CD",
    "Sm"   = "#0000EE",

    "ALe5"   = "#8B5742",
    "B51"    = "#CD8162",
    "B10to50"= "#EE9572",
    "Old50"  = "#FFA07A",

    "Hi45" = "#8B8B00",
    "Lo"   = "#CDCD00",
    "No"   = "#FFFF00",

    "HaVe" = "#8B475D",
    "Stro" = "#CD6889",
    "Npro" = "#EE799F",
    "Fra"  = "#DB7093"

  )
}

#' Color Palette for Freshwater Functional Traits
#'
#' Returns a named vector with predefined colors for functional trait modalities
#' of freshwater benthic macroinvertebrates.
#'
#' @return A named vector with color codes for functional traits including feeding
#'   habits, microhabitat preferences, substrate type, food sources, feeding type,
#'   locomotion, respiration, body size, life cycle duration, and voltinism.
#'
#' @details
#' This comprehensive palette covers multiple functional trait categories:
#' \itemize{
#'   \item \strong{Feeding habits (FED):} gra (grazers), min (miners), xyl (xylophagous),
#'     shr (shredders), gat (gatherers), aff (active filter feeders), pff (passive filter feeders),
#'     pre (predators), par (parasites), othF (other)
#'   \item \strong{Microhabitat (MH):} arg (argyllal), pel (pelal), psa (psammal),
#'     aka (akal), lit (lithal), phy (phytal), pom (POM), othM (other)
#'   \item \strong{Substrate (SBT):} fbcp (flags/boulders/cobbles/pebbles), grvl (gravel),
#'     sand, silt, macp (macrophytes), micp (microphytes), twro (twigs/roots),
#'     odli (organic detritus), mud, othS (other)
#'   \item \strong{Food (FD):} mior (microorganisms), detl1 (detritus <1mm), dpg1 (dead plants),
#'     limi (living microphytes), lima (living macrophytes), dag1 (dead animals),
#'     lmic (living microinvertebrates), lmac (living macroinvertebrates), vert (vertebrates)
#'   \item \strong{Feeding type (FEH):} abs (absorbers), dpf (deposit feeders), shr (shredders),
#'     scr (scrapers), fif (filter-feeders), pie (piercers), preFt (predators), parFt (parasites)
#'   \item \strong{Locomotion type (LoTy):} sws (swimming/skating), swd (swimming/diving),
#'     bub (burrowing/boring), spw (sprawling/walking), ses (sessile), othL (other)
#'   \item \strong{Locomotion substrate (LoSb):} flr (flying), ssw (surface swimming),
#'     wsw (water swimming), crw (crawling), bur (burrowing), int (interstitial),
#'     tat (temporarily attached), pat (permanently attached)
#'   \item \strong{Respiration (Res):} teg (tegument), gil (gill), pls (plastron),
#'     spi (spiracle), ves (hydrostatic vesicle), sur (surface excursion)
#'   \item \strong{Body size (MPS):} L25 (≤0.25cm), 25To5 (>0.25-0.5cm), 5To1 (>0.5-1cm),
#'     1To2 (>1-2cm), 2To4 (>2-4cm), 4To8 (>4-8cm), H8 (>8cm)
#'   \item \strong{Life cycle (LiCy):} <1y (≤1 year), >1y (>1 year)
#'   \item \strong{Voltinism (Vol):} <1 (semivoltine), 1 (monovoltine), >1 (polyvoltine)
#' }
#'
#' Color schemes are organized by functional categories using distinct color gradients
#' to facilitate visual interpretation of community functional structure.
#'
#' @examples
#' # Get the complete freshwater trait palette
#' fw_colors <- FFCol()
#'
#' # View all available trait modalities
#' names(fw_colors)
#'
#' # Get colors for feeding habits only
#' feeding_colors <- fw_colors[c("gra", "min", "xyl", "shr", "gat",
#'                                "aff", "pff", "pre", "par")]
#'
#' # Use in a plot
#' library(ggplot2)
#' trait_data <- data.frame(
#'   trait = c("gra", "shr", "pre"),
#'   abundance = c(35, 28, 15)
#' )
#' ggplot(trait_data, aes(x = trait, y = abundance, fill = trait)) +
#'   geom_bar(stat = "identity") +
#'   scale_fill_manual(values = FFCol())
#'
#' @export
FFCol <- function() {
  c(
    "gra" = "#00441B",
    "min" = "#006D2C",
    "xyl" = "#238B45",
    "shr" = "#41AB5D",
    "gat" = "#74C476",
    "aff" = "#A1D99B",
    "pff" = "#C7E9C0",
    "pre" = "#31A354",
    "par" = "#2C7FB8",
    "othF" = "#EDF8E9",
    "arg" = "#08306B",
    "pel" = "#08519C",
    "psa" = "#2171B5",
    "aka" = "#4292C6",
    "lit" = "#6BAED6",
    "phy" = "#9ECAE1",
    "pom" = "#C6DBEF",
    "othM" = "#DEEBF7",
    "fbcp" = "#3E2723",
    "grvl" = "#5D4037",
    "sand" = "#795548",
    "silt" = "#8D6E63",
    "macp" = "#A1887F",
    "micp" = "#BCAAA4",
    "twro" = "#6D4C41",
    "odli" = "#4E342E",
    "mud"  = "#3E2723",
    "othS" = "#D7CCC8",
    "mior" = "#7F2704",
    "detl1" = "#A63603",
    "dpg1" = "#D94801",
    "limi" = "#F16913",
    "lima" = "#FD8D3C",
    "dag1" = "#FDAE6B",
    "lmic" = "#FDD0A2",
    "lmac" = "#FEE6CE",
    "vert" = "#E6550D",
    "abs" = "#67000D",
    "dpf" = "#A50F15",
    "shr" = "#CB181D",
    "scr" = "#EF3B2C",
    "fif" = "#FB6A4A",
    "pie" = "#FC9272",
    "preFt" = "#DE2D26",
    "parFt" = "#A50F15",
    "sws" = "#3F007D",
    "swd" = "#54278F",
    "bub" = "#6A51A3",
    "spw" = "#807DBA",
    "ses" = "#9E9AC8",
    "othL" = "#DADAEB",
    "flr" = "#00441B",
    "ssw" = "#006D2C",
    "wsw" = "#238B45",
    "crw" = "#41AB5D",
    "bur" = "#005A32",
    "int" = "#238B45",
    "tat" = "#66C2A4",
    "pat" = "#99D8C9",
    "teg" = "#0868AC",
    "gil" = "#2B8CBE",
    "pls" = "#4EB3D3",
    "spi" = "#7BCCC4",
    "ves" = "#A8DDB5",
    "sur" = "#CCEBC5",
    "L25" = "#FFFFCC",
    "25To5" = "#FFEDA0",
    "5To1" = "#FED976",
    "1To2" = "#FEB24C",
    "2To4" = "#FD8D3C",
    "4To8" = "#F03B20",
    "H8" = "#BD0026",
    "<1y" = "#636363",
    ">1y" = "#BDBDBD",
    "<1" = "#08589E",
    "1" = "#2B8CBE",
    ">1" = "#7BCCC4"
  )
}
