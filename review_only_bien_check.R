library(BIEN)
library(tidyverse)

# start with species list
species.list <- read.csv("data/species_list_bien.csv") %>% 
  dplyr::select(SpeciesBinomial)

species.list

species.traits <- BIEN_trait_species(species=c(
  "Acer circinatum",
  "Acer glabrum",
  "Achillea millefolium",
  "Amelanchier alnifolia",
  "Arctostaphylos uva-ursi",
  "Athyrium filix-femina",
  "Cassiope mertensiana",
  "Calamagrostis rubescens",
  "Ceanothus velutinus",
  "Chamaenerion angustifolium",
  "Chimaphila umbellata",
  "Clintonia uniflora",
  "Cornus canadensis",
  "Cornus sericea",
  "Gaultheria ovatifolia",
  "Gaultheria shallon",
  "Goodyera oblongifolia",
  "Gymnocarpium dryopteris",
  "Hieracium albiflorum",
  "Holodiscus discolor",
  "Linnaea borealis",
  "Luetkea pectinata",
  "Mahonia nervosa",
  "Menziesia ferruginea",
  "Oplopanax horridus",
  "Paxistima myrsinites",
  "Phyllodoce empetriformis",
  "Polystichum munitum",
  "Pteridium aquilinum",
  "Rhododendron albiflorum",
  "Rubus lasiococcus",
  "Rubus parviflorus",
  "Rubus pedatus",
  "Rubus spectabilis",
  "Sorbus sitchensis",
  "Spiraea betulifolia",
  "Tiarella trifoliata",
  "Trientalis borealis",
  "Vaccinium alaskaense",
  "Vaccinium deliciosum",
  "Vaccinium membranaceum",
  "Valeriana sitchensis"))

traits.available <- unique(species.traits$trait_name)
