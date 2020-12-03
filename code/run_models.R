#devtools::install_github("pbs-assess/sdmTMB")
library(ggplot2)
library(raster)
library(dplyr)
library(sdmTMB)
library(sp)

# haul data includes environmental covariates with location information
haul = nwfscSurvey::PullHaul.fn(YearRange = c(2003, 2018), SurveyName = "NWFSC.Combo")

# read in the grid cell data from the survey design
grid_cells = readxl::read_excel("data/Selection Set 2018 with Cell Corners.xlsx")
grid_cells = dplyr::mutate(grid_cells,
  depth_min = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[1]),
  depth_max = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[2]))

# Figure out the grid cell corresponding to each tow location
haul$Cent.Lat = NA
haul$Cent.Lon = NA
haul$Cent.ID = NA
for(i in 1:nrow(haul)) {
  indx = which(grid_cells$NW.LAT > haul$latitude_dd[i] &
      grid_cells$SW.LAT < haul$latitude_dd[i] &
      grid_cells$NW.LON < haul$longitude_dd[i] &
      grid_cells$NE.LON > haul$longitude_dd[i])
  if(length(indx) > 0) {
    haul$Cent.ID[i] = grid_cells$Cent.ID[indx]
    haul$Cent.Lat[i] = grid_cells$Cent.Lat[indx]
    haul$Cent.Lon[i] = grid_cells$Cent.Long[indx]
  }
}

# project lat/lon to UTM, after removing missing values and unsatisfactory hauls
haul = haul %>% filter(!is.na(Cent.Lon), performance == "Satisfactory")

haul_trans = haul
coordinates(haul_trans) <- c("Cent.Lon", "Cent.Lat")
proj4string(haul_trans) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
haul_trans <- spTransform(haul_trans, CRS(newproj))
haul_trans = as.data.frame(haul_trans)
haul_trans$Cent.Lon = haul_trans$Cent.Lon/10000
haul_trans$Cent.Lat = haul_trans$Cent.Lat/10000
haul_trans$year = as.numeric(substr(haul_trans$date_yyyymmdd,1,4))

haul$X = haul_trans$Cent.Lon
haul$Y = haul_trans$Cent.Lat
haul$year = haul_trans$year

# center and scale depth, removing NAs
haul = dplyr::filter(haul, !is.na(depth_hi_prec_m))
haul$log_depth_scaled = scale(log(haul$depth_hi_prec_m))
haul$log_depth_scaled2 = haul$log_depth_scaled ^ 2

coordinates(grid_cells) <- c("Cent.Long", "Cent.Lat")
proj4string(grid_cells) <- CRS("+proj=longlat +datum=WGS84")
grid_cells <- spTransform(grid_cells, CRS(newproj))

# make prediction raster roughly from grid_cell centroids, given standard cell dimensions (here in meters, converted from nm)
predict_raster = raster(grid_cells, resolution = c(2778,3704), vals = NULL)
## load custom bathymetry raster
bathy_hiRes <- raster("data/bathy_clipped")
bathy_hiRes <- bathy_hiRes / 10 # units were originally decimeters, so convert to meters
# aggregate and project bathymetry to survey grid cells, the absolute minimum resolution of the prediction grid
bathy_raster <- projectRaster(bathy_hiRes, predict_raster, crs = newproj, method="bilinear")
# load Cowcod Conservation Areas, not included in trawl survey, and reproject
CCA = rgdal::readOGR('data/kv299cy7357.shp')
CCA = sp::spTransform(CCA, sp::CRS(newproj))
# mask CCA from bathymetry raster used for prediction
bathy_raster = suppressWarnings(raster::mask(bathy_raster, CCA, inverse = TRUE))
# create matrix of point data with coordinates and depth from raster
wc_grid <- as.data.frame(rasterToPoints(bathy_raster))
colnames(wc_grid) = c("X", "Y", "depth")

# scale covariates
wc_grid$log_depth_scaled <- (log(wc_grid$depth * -1) - mean(log(haul$depth_hi_prec_m))) / sd(log(haul$depth_hi_prec_m))
wc_grid$log_depth_scaled2 <- wc_grid$log_depth_scaled ^ 2
wc_grid$X <- wc_grid$X/10000
wc_grid$Y <- wc_grid$Y/10000

saveRDS(wc_grid, file=paste0("data/wc_grid.rds")) # save prediction grid


# catch data includes catch, effort, etc. This takes a few minutes to grab all ~ 900 spp
catch = nwfscSurvey::PullCatch.fn(YearRange = c(2003, 2018), SurveyName="NWFSC.Combo")
# format to later join catch and haul
names(catch) = tolower(names(catch))
catch$trawl_id = as.numeric(catch$trawl_id)

catch$common_name = NA
catch$common_name[which(catch$scientific_name=="Microstomus pacificus")] = "Dover sole"
catch$common_name[which(catch$scientific_name=="Sebastolobus alascanus")] = "shortspine thornyhead"
catch$common_name[which(catch$scientific_name=="Sebastolobus altivelis")] = "longspine thornyhead"
catch$common_name[which(catch$scientific_name=="Atheresthes stomias")] = "arrowtooth flounder"
catch$common_name[which(catch$scientific_name=="Eopsetta jordani")] = "petrale sole"
catch$common_name[which(catch$scientific_name=="Ophiodon elongatus")] = "lingcod"
catch$common_name[which(catch$scientific_name=="Hippoglossus stenolepis")] = "Pacific halibut" # few obs
catch$common_name[which(catch$scientific_name=="Glyptocephalus zachirus")] = "rex sole"
catch$common_name[which(catch$scientific_name=="Parophrys vetulus")] = "English sole"
catch$common_name[which(catch$scientific_name=="Anoplopoma fimbria")] = "sablefish"
catch$common_name[which(catch$scientific_name=="Sebastes pinniger")] = "canary rockfish" # few obs
catch$common_name[which(catch$scientific_name=="Sebastes flavidus")] = "yellowtail rockfish" # few obs
catch$common_name[which(catch$scientific_name=="Sebastes entomelas")] = "widow rockfish"
catch$common_name[which(catch$scientific_name=="Sebastes ruberrimus")] = "yelloweye rockfish" # few obs
catch$common_name[which(catch$scientific_name=="Sebastes diploproa")] = "splitnose rockfish"
catch$common_name[which(catch$scientific_name=="Sebastes alutus")] = "Pacific ocean perch"
catch$common_name[which(catch$scientific_name=="Sebastes crameri")] = "darkblotched rockfish"
catch$common_name[which(catch$scientific_name=="Sebastes paucispinis")] = "bocaccio" # few obs
catch$common_name[which(catch$scientific_name=="Sebastes goodei")] = "chilipepper" # few obs
catch$common_name[which(catch$scientific_name=="Gadus macrocephalus")] = "Pacific cod" # few obs
catch$common_name[which(catch$scientific_name=="Raja rhina")] = "longnose skate"
catch$common_name[which(catch$scientific_name=="Raja binoculata")] = "big skate"
catch$common_name[which(catch$scientific_name=="Squalus suckleyi")] = "spiny dogfish"
catch$common_name[which(catch$scientific_name=="Citharichthys sordidus")] = "Pacific sanddab"
catch$common_name[which(catch$scientific_name=="Hydrolagus colliei")] = "spotted ratfish"

catch = dplyr::filter(catch, !is.na(common_name))

# Loop over species
species = unique(catch$common_name)

use_cv = FALSE

for(spp in 1:length(species)) {
  subset = dplyr::filter(catch, common_name == species[spp])
  haul_new = dplyr::left_join(haul, subset)

  # iterate fits over a range of number of knots,
  # using Tweedie predictive density to evaluate performance
  performance <- data.frame(
    knots = seq(350,450, by = 100),
    tweedie_dens = 0
  )

  for (k in seq_len(nrow(performance))) {

    # when using cross-validation, create fold ids based on latitude quantiles
    haul_new$fold = 5
    haul_new$fold[which(haul_new$lat < quantile(haul_new$lat,0.8))] = 4
    haul_new$fold[which(haul_new$lat < quantile(haul_new$lat,0.6))] = 3
    haul_new$fold[which(haul_new$lat < quantile(haul_new$lat,0.4))] = 2
    haul_new$fold[which(haul_new$lat < quantile(haul_new$lat,0.2))] = 1

    # Set NA CPUEs to 0
    haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] = 0

    if(use_cv==TRUE) {
    density_model <- sdmTMB_cv(formula = cpue_kg_km2 ~ log_depth_scaled + log_depth_scaled2 + as.factor(year),
      data = haul_new, x = "X", y = "Y", k_folds=max(haul_new$fold), fold_ids = "fold",
      n_knots = performance$knots[k],
      time = "year", anisotropy = TRUE,
      silent = TRUE, spatial_trend = TRUE, spatial_only = FALSE,
      family = tweedie(link = "log"),
      control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

    if(!dir.exists(paste0("results/",species[spp]))) dir.create(paste0("results/",species[spp]))
    saveRDS(density_model, file=paste0("results/",species[spp],"/",species[spp],"_",performance$knots[k],"_density_all_yr.rds"))
    } else {
      c_spde <- make_spde(haul_new$X, haul_new$Y, n_knots = performance$knots[k])
      density_model <- sdmTMB(formula = cpue_kg_km2 ~ log_depth_scaled + log_depth_scaled2 + as.factor(year),
        data = haul_new,
        time = "year", spde = c_spde, anisotropy = TRUE,
        silent = TRUE, spatial_trend = TRUE, spatial_only = FALSE,
        family = tweedie(link = "log"),
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

      if(!dir.exists(paste0("results/",species[spp]))) dir.create(paste0("results/",species[spp]))
      saveRDS(density_model, file=paste0("results/",species[spp],"/",species[spp],"_",performance$knots[k],"_density_all_yr.rds"))

    }

    # validate against the test set
    performance[k, "tweedie_dens"] = sum(density_model$data$cv_loglik)
    saveRDS(performance, file = paste0("results/",species[spp],"_performance_yr.rds"))
  }
}
