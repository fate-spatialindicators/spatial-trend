library(sdmTMB)
library(ggplot2)
library(viridis)
library(dplyr)
library(fpc)
library(gridExtra)

# choose which models based on # of knots
knots = 450

# select years to predict and gather replicate prediction data for each of these years
years = 2003:2018
wc_grid = readRDS("data/wc_grid.rds")
wc_grid_all = wc_grid
wc_grid_all$year = 2003

for(i in 2:length(unique(years))) {
  wc_grid$year = unique(years)[i]
  wc_grid_all = rbind(wc_grid_all, wc_grid)
}

species = dir("results")[grep("performance_yr.rds", dir("results"))]
species = unlist(lapply(strsplit(species,"_"), getElement, 1))
species_drop = c("chilipepper", "longspine thornyhead", "Pacific cod", "Pacific sanddab",
                 "yelloweye rockfish", "yellowtail rockfish")# species with poor model convergence/fit
species = species[!(species %in% species_drop)]

trend_cluster_plots = list()
intercept_cluster_plots = list()
all_clust = NULL

# loop over species
for(spp in 1:length(species)) {
     d = readRDS(paste0("results/", species[spp],
       "/",species[spp],"_",
       knots,"_density_all_yr.rds"))
     p = predict(d, newdata=wc_grid_all)

     # apply post-hoc clustering based on the trend/intercept and latitude, for 1 year (trend/ intercept same for all years)
     unique_lon = unique(p$X)[seq(1,length(unique(p$X)),by=10)] # downsampling longitude to speed things up
     # We use the scale function to standardize lat/trend, otherwise lat totally dominates
     trend_clust = pamk(scale(dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("Y","zeta_s")]))
     d2003 = dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("X","Y","zeta_s")]
     d2003$trend_cluster = trend_clust$pamobject$clustering

     # compile data for stripplot of trend clusters by latitude, colored by trend anomaly, for all species
     d2003$zeta_s = scale(d2003$zeta_s, scale = FALSE) # mean center slopes to remove influence of total population abundance trend
     all_clust = rbind(all_clust,
                           d2003 %>% group_by(trend_cluster) %>% mutate(mean_zeta_s = mean(zeta_s), species = species[spp]))

     #d2003$trend_cluster = as.factor(d2003$trend_cluster)
}

# save results
save(all_clust, file = "results/plotdata_stripplot.Rdata")

# make stripplot of clusters by latitude, colored by mean slope
ggplot(all_clust, aes(x=species, y=Y, group = trend_cluster, color=cut(mean_zeta_s, c(-Inf, -0.01, 0.01, Inf)))) +
    geom_jitter(position=position_dodge(0.4), size=0.5) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("Northings (10s km)") + geom_hline(yintercept=c(450,385)) +
    scale_color_manual(name = "trend anomaly",
                       values = c("(-Inf,-0.01]" = "red",
                                  "(-0.01,0.01]" = "darkgrey",
                                  "(0.01, Inf]" = "blue"),
                       labels = c("-", "0", "+"))
ggsave(filename = "plots/Fig4_trendcluster_stripplot_450.pdf",
       width = 10, height = 6, units = c("in"))
