devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
#devtools::install_github("pbs-assess/sdmTMB")
library(sdmTMB)
library(ggplot2)
library(viridis)
library(dplyr)
library(fpc)
library(gridExtra)
library(rnaturalearth)

# choose which models based on # of knots
knots = 350
# choose year to use as reference example (in figure 3)
ex_year = 2011

# select years to predict and gather replicate prediction data for each of these years
years = 2003:2018
wc_grid = readRDS("data/wc_grid.rds")
wc_grid_all = wc_grid
wc_grid_all$year = 2003

for(i in 2:length(unique(years))) {
  wc_grid$year = unique(years)[i]
  wc_grid_all = rbind(wc_grid_all, wc_grid)
}

# create subsets of grid for biogeographic regions
wc_grid_all_N = filter(wc_grid_all, Y > 450)
wc_grid_all_C = filter(wc_grid_all, Y <= 450, Y >= 385)
wc_grid_all_S = filter(wc_grid_all, Y < 385)

# define species to plot
species = dir("results")[grep("performance_yr.rds", dir("results"))]
species = unlist(lapply(strsplit(species,"_"), getElement, 1))
species_drop = c("chilipepper", "longspine thornyhead", "Pacific cod", "Pacific sanddab",
                 "yelloweye rockfish", "yellowtail rockfish")# species with poor model convergence/fit
species = species[!(species %in% species_drop)]

anisotropy_plots = list()
qq_plots = list()
residuals_plots = list()
prediction_plots = list()
prediction_plots_ex = list()
spatiotemporal_plots = list()
trend_plots = list()
mean_dens_plots = list()
cluster_plots = list()
COG_plots = list()
all_clust = NULL

# plotting functions
plot_map_point <- function(dat, column = "omega_s") {
  ggplot(dat, aes_string("X", "Y", colour = column)) +
    geom_point(size=0.1) +
    xlab("Eastings") +
    ylab("Northings")
}
plot_map_raster <- function(dat, column = "omega_s") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    xlab("Eastings") +
    ylab("Northings")
}


# loop over species
for(spp in 1:length(species)) {
  d = readRDS(paste0("results/", species[spp],
                     "/",species[spp],"_",
                     knots,"_density_all_yr.rds"))
  # below 2 lines necessary for models fit with older sdmTMB versions
  d$tmb_data$weights_i = rep(1, length(d$tmb_data$y_i))
  d$tmb_data$calc_quadratic_range = as.integer(FALSE)

  p = predict(d, newdata=wc_grid_all)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # apply post-hoc clustering based on the trend and latitude, for 1 year (trend/ intercept same for all years)

  # Original method, but including testing whether 1 cluster is supported
  #unique_lon = unique(p$X)[seq(1,length(unique(p$X)),by=10)] # downsampling longitude to speed things up
  #clust = pamk(scale(dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("Y","zeta_s")]), krange = 1:10)
  #d2003 = dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("X","Y","zeta_s")]

  # Alternative method 1: as in main manuscript figure but without latitude
  #unique_lon = unique(p$X)[seq(1,length(unique(p$X)),by=10)] # downsampling longitude to speed things up
  #clust = pamk(dplyr::filter(p, year == 2003, X %in% unique_lon)[,"zeta_s"])
  #d2003 = dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("X","Y","zeta_s")]

  # Alternative method 1b: as in main manuscript figure but without latitude, constrained to 2-3 clusters
  unique_lon = unique(p$X)[seq(1,length(unique(p$X)),by=10)] # downsampling longitude to speed things up
  clust = pamk(dplyr::filter(p, year == 2003, X %in% unique_lon)[,"zeta_s"], krange = 2:3)
  d2003 = dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("X","Y","zeta_s")]

  # Alternative method 2: as method 1 but without downscaling
  #clust = pamk(dplyr::filter(p, year == 2003)[,"zeta_s"], usepam = FALSE)
  #d2003 = dplyr::filter(p, year == 2003)[,c("X","Y","zeta_s")]

  d2003$cluster = clust$pamobject$clustering

  # compile data for stripplot of trend clusters by latitude, colored by trend anomaly, for all species
  d2003$zeta_s = scale(d2003$zeta_s, scale = FALSE) # mean center slopes to remove influence of total population abundance trend
  all_clust = rbind(all_clust,
                    d2003 %>% group_by(cluster) %>% mutate(mean_zeta_s = mean(zeta_s), species = species[spp]))

  d2003$cluster = as.factor(d2003$cluster)

  # cluster maps
  cluster_plots[[spp]] = plot_map_point(d2003, "cluster") +
    theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
          axis.text = element_blank(), axis.title = element_blank(),
          legend.position = "none") +
    labs(title = "Trend cluster")

}

# save results
#save.image(file = paste0("results/plotdata_trendcluster_downscaled_",knots,".Rdata"))
#save.image(file = paste0("results/plotdata_trendonlycluster_downscaled_",knots,".Rdata"))
save.image(file = paste0("results/plotdata_trendonlycluster_downscaled_clust_lim_",knots,".Rdata"))
#save.image(file = paste0("results/plotdata_trendonlycluster_",knots,".Rdata"))


# make stripplot of clusters by latitude, colored by mean slope
ggplot(all_clust, aes(x=species, y=Y, group = cluster, color=cut(mean_zeta_s, c(-Inf, -0.01, 0.01, Inf)))) +
  geom_jitter(position=position_dodge(0.4), size=0.5) +
  theme_sleek()+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_text(size=13),
        legend.text = element_text(size=14)) +
  ylab("Northings (10s km)") +
  geom_hline(yintercept=c(450,385)) +
  scale_color_manual(name = "Trend anomaly",
                     values = c("(-Inf,-0.01]" = "red",
                                "(-0.01,0.01]" = "darkgrey",
                                "(0.01, Inf]" = "blue"),
                     labels = c("-", "0", "+")) +
  guides(color = guide_legend(override.aes = list(size=5)))
ggsave(filename = paste0("plots/Fig4_trendonlycluster_stripplot_clust_lim",knots,"_0.01.pdf"),
       width = 10, height = 6, units = c("in"))


# analyze how likely breaks are to be close to biogeographic breaks
clust_breaks <- ungroup(all_clust) %>%
  group_by(species, cluster) %>%
  mutate(low_break = min(Y), high_break = max(Y)) %>%
  ungroup() %>%
  distinct(species, cluster, .keep_all = TRUE) %>%
  group_by(species) %>%
  mutate(match_breaks = low_break <= 395 & low_break >= 375 | low_break <= 460 & low_break >= 440 |
           high_break <= 395 & high_break >= 375 | high_break <= 460 & high_break >= 440,
         between_breaks = low_break > 395 & low_break < 440 | high_break > 395 & high_break < 440)

clust_breaks %>% summarize(prop_match = mean(match_breaks), prop_between = mean(between_breaks)) # proportions within species
ungroup(clust_breaks) %>% summarize(prop_match = mean(match_breaks), prop_between = mean(between_breaks)) # overall proportions



