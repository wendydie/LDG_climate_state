# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: options.R
# Last updated: 2026-08-01
# -----------------------------------------------------------------------

# Parameters and options for analyses
params <- list(
  # Download updated dataset
  download = FALSE,
  # PBDB base url
  base_url = "https://paleobiodb.org/data1.2/occs/list.csv",
  # PBDB Download options
  query = list(
    base_name = "Animalia^Vertebrata",
    taxon_reso = "genus",
    idtype = "latest",
    taxon_status = "valid",
    idqual = "genus_certain",
    pres = "regular",
    interval = "Tremadocian,Holocene",
    envtype = "marine",
    show = "genus,pres,strat,coll,coords,loc,class,env,resgroup"),
  # The geological rank for conducting analyses
  rank = "stage",
  # The Geological Time Scale to be used
  GTS = "international ages",
  # Naming convention for temporal bin midpoint
  age = "bin_midpoint",
  # How should occurrences be temporally binned?
  method = "majority",
  # Threshold for majority binning rule
  threshold = 50,
  # Naming convention for longitude
  lng = "lng",
  # Naming convention for latitude
  lat = "lat",
  # Naming convention for palaeolongitude
  p_lng = "p_lng",
  # Naming convention for palaeolatitude
  p_lat = "p_lat",
  # Number of lat/lng decimal places to define stacked collections
  n_decs = 2,
  # Which Global Plate Models should be used?
  models = c("PALEOMAP"),
  # Which spacing should be used for spatial binning (in km)?
  spacing = 250,
  # Whether to clean the data again depends on the situation.
  clean_again = FALSE,
  # Coordinate reference system, default is WGS84
  crs = 'epsg:4326', 
  # Diversity order(s) to calculate; 0 for species richness
  q = c(0),
  # Data type; 'incidence_freq' for presence-absence data aggregated into frequencies
  datatype = "incidence_freq", 
  # Standardization base; 'coverage' ensures comparisons at the same coverage level
  base = "coverage", 
  # Desired sample coverage level (e.g., 0.7 = 70% of species expected to be detected)
  level = 0.7, 
  # Number of bootstrap replicates for confidence intervals; 0 means no bootstrapping
  nboot = 50,
  notify = TRUE
)


