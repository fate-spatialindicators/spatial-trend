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

     # make plot of predictions from full model (all fixed + random effects)
     prediction_plots[[spp]] = plot_map_raster(p, "est") +
       facet_wrap(~year) +
       coord_fixed() +
       scale_fill_viridis_c() +
       ggtitle(paste0(species[spp],"_predicted_density"))
     # make plot of predictions from spatiotemporal random effects
     spatiotemporal_plots[[spp]] = plot_map_raster(p, "epsilon_st") +
       facet_wrap(~year) +
       coord_fixed() +
       scale_fill_viridis_c() +
       ggtitle(paste0(species[spp],"_ST"))

     # make plot of general spatial pattern of theta
     if(spp != length(species)){
       trend_plots[[spp]] = plot_map_raster(dplyr::filter(p,year==2003), "zeta_s") +
         scale_fill_gradient2() +
         theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
             axis.title.x = element_blank(),
             axis.text = element_blank(),
             legend.key.width = unit(0.1,"cm"),
             legend.title = element_blank(),
             legend.position = "right",
             legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
         labs(title = "Trend", y = species[spp])
     }else{
       trend_plots[[spp]] = plot_map_raster(dplyr::filter(p,year==2003), "zeta_s") +
         scale_fill_gradient2() +
         theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
               axis.title.x = element_blank(),
               legend.key.width = unit(0.1,"cm"),
               legend.title = element_blank(),
               legend.position = "right",
               legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
         labs(title = "Trend", y = species[spp])
     }

     # make plot of predictions from full model for example year
     p$est_exp <- exp(p$est)
     prediction_plots_ex[[spp]] = plot_map_raster(dplyr::filter(p,year==ex_year), "est_exp") +
       scale_fill_viridis_c() +
       theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
             axis.title = element_blank(),
             axis.text = element_blank(),
             legend.key.width = unit(0.2,"cm"),
             legend.title = element_blank(),
             legend.position = "right",
             legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
       labs(title = "Density 2011")

     ## plot mean density
     # extract all coefficients
     sr_se <- summary(d$sd_report)[,"Std. Error"]
     b_j <- unname(d$model$par[grep("b_j", names(d$model$par))])
     b_j_se <- unname(sr_se[grep("b_j", names(sr_se))])
     mm <- cbind(b_j, b_j_se)
     colnames(mm) <- c("coef.est", "coef.se")
     row.names(mm) <- colnames(model.matrix(d$formula, d$data))
     # get mean density
     p$mean_exp <- exp(dplyr::filter(p,year==2003)$omega_s + 0.5*dplyr::filter(p,year==2003)$zeta_s + b_j[1] + b_j[10] + p$log_depth_scaled*b_j[2] + p$log_depth_scaled2*b_j[3])

     mean_dens_plots[[spp]] = plot_map_raster(dplyr::filter(p,year==ex_year), "mean_exp") +
       scale_fill_viridis_c() +
       theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
             axis.title = element_blank(),
             axis.text = element_blank(),
             legend.key.width = unit(0.2,"cm"),
             legend.title = element_blank(),
             legend.position = "right",
             legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
       labs(title = "Mean density")

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # apply post-hoc clustering based on the trend and latitude, for 1 year (trend/ intercept same for all years)
     unique_lon = unique(p$X)[seq(1,length(unique(p$X)),by=10)] # downsampling longitude to speed things up
     # We use the scale function to standardize lat/trend, otherwise lat totally dominates
     clust = pamk(scale(dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("Y","zeta_s")]))
     d2003 = dplyr::filter(p, year == 2003, X %in% unique_lon)[,c("X","Y","zeta_s")]
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

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # make timeseries plots of COG of latitude, contrasting coastwide vs for each biogeographic region
     p_All = predict(d, newdata=wc_grid_all, return_tmb_object = TRUE)
     p_N = predict(d, newdata=wc_grid_all_N, return_tmb_object = TRUE)
     p_C = predict(d, newdata=wc_grid_all_C, return_tmb_object = TRUE)
     p_S = predict(d, newdata=wc_grid_all_S, return_tmb_object = TRUE)
     COG_All = get_cog(p_All)
     COG_N = get_cog(p_N) %>% mutate(region = "N")
     COG_C = get_cog(p_C) %>% mutate(region = "C")
     COG_S = get_cog(p_S) %>% mutate(region = "S")
     # get proportion of density for each region
     sum_D = sum(exp(p$est)) # sum of density coastwide
     COG_N$prop_D = sum(exp(p_N$data$est)) / sum_D
     COG_C$prop_D = sum(exp(p_C$data$est)) / sum_D
     COG_S$prop_D = sum(exp(p_S$data$est)) / sum_D

     # set up structures for plotting COG in dimension indicated by axis, where prop_D is proportion of density in each region
     axis = "Y"
     COG = COG_All %>% filter(coord == axis) %>% mutate(prop_D = 1) # coastwide
     COG_region = rbind(COG_N,COG_C,COG_S) %>% filter(coord == axis) # by region

     # add x axis labels and legend only to specific panels to be combined later
     if(!(spp %in% c(17,19))){
       COG_plots[[spp]] = ggplot(COG_region, aes(year, est)) +
         geom_ribbon(data = COG, aes(ymin = lwr, ymax = upr), fill = "grey70") +
         geom_line(data = COG, color = "black") +
         geom_ribbon(aes(group = region, fill = prop_D, ymin = lwr, ymax = upr), alpha = 0.5) +
         geom_line(aes(group = region, color = prop_D)) +
         scale_color_gradient(low = "#feb24c", high = "#a50f15", breaks = c(0,1), limits = c(0, 1)) +
         scale_fill_gradient(low = "#feb24c", high = "#a50f15", breaks = c(0,1), limits = c(0, 1)) +
         labs(title = species[spp], tag = LETTERS[spp]) +
         theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
               axis.text.x = element_blank(), plot.margin = unit(c(0,0,1,3), "pt"),
               legend.position = "none")
     }
      else{
        COG_plots[[spp]] = ggplot(COG_region, aes(year, est)) +
          geom_ribbon(data = COG, aes(ymin = lwr, ymax = upr), fill = "grey70") +
          geom_line(data = COG, color = "black") +
          geom_ribbon(aes(group = region, fill = prop_D, ymin = lwr, ymax = upr), alpha = 0.5) +
          geom_line(aes(group = region, color = prop_D)) +
          scale_color_gradient(low = "#feb24c", high = "#a50f15", breaks = c(0,1), limits = c(0, 1)) +
          scale_fill_gradient(low = "#feb24c", high = "#a50f15", breaks = c(0,1), limits = c(0, 1)) +
          labs(title = species[spp], tag = LETTERS[spp]) +
          theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
                plot.margin = unit(c(0,0,0,3), "pt"),
                legend.position = "none")
      }
}

# save results
#save.image(file = paste0("results/plotdata_lim_yearfixed_",knots,".Rdata"))

# COG timeseries plots
ggsave(filename = paste0("plots/FigS2_COG_color_",knots,".pdf"),
       plot = arrangeGrob(grobs = COG_plots, ncol = 5, bottom = "year",
                          left = grid::textGrob("COG Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 12, height = 8, units = c("in"))

# plot of predictions from full model (all fixed + random effects)
ggsave(filename = paste0("plots/predicted_density_maps_",knots,".pdf"),
       plot = marrangeGrob(prediction_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))

# plot only spatiotemporal random effects
ggsave(filename = paste0("plots/st_maps_",knots,".pdf"),
       plot = marrangeGrob(spatiotemporal_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))

# make stripplot of clusters by latitude, colored by mean slope
ggplot(all_clust, aes(x=species, y=Y, group = cluster, color=cut(mean_zeta_s, c(-Inf, -0.01, 0.01, Inf)))) +
  geom_jitter(position=position_dodge(0.4), size=0.5) +
  theme_sleek()+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Northings (10s km)") + geom_hline(yintercept=c(450,385)) +
  scale_color_manual(name = "Trend anomaly",
                     values = c("(-Inf,-0.01]" = "red",
                                "(-0.01,0.01]" = "darkgrey",
                                "(0.01, Inf]" = "blue"),
                     labels = c("-", "0", "+"))
ggsave(filename = paste0("plots/Fig4_trendcluster_stripplot_",knots,"_0.01.pdf"),
       width = 10, height = 6, units = c("in"))


## plot COG along with trend, trend cluster, and density, for selected species subsets
spp_subset1 = c(1, 3, 7, 12, 14, 15)
spp_subset2 = (1:spp)[!(1:spp %in% spp_subset1)]
trend_plots1 = list()
mean_dens_plots1 = list()
cluster_plots1 = list()
prediction_plots_ex1 = list()
COG_plots1 = list()

# get coastlines and transform to projection and scale of data
shore <- rnaturalearth::ne_countries(continent = "north america", scale = "medium", returnclass = "sp")
shore <- sp::spTransform(shore, CRS = "+proj=utm +zone=10 ellps=WGS84")
shore <- fortify(shore)
shore$long <- shore$long/10000
shore$lat <- shore$lat/10000

# use customized theme
theme_map <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      #panel.grid       = element_line(colour = "grey87"),
      #panel.grid.major = element_line(size = rel(0.5)),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(fill = NA, colour = "grey70", size = rel(1)),
      axis.ticks       = element_line(colour = "grey70", size = rel(0.5)),
    )
}
theme_set(theme_map())

for(spp in 1:length(species)){
  trend_plots1[[spp]] = trend_plots[[spp]] + annotation_map(shore, color = "black", fill = "white", size=0.1) +
    theme(text = element_text(size = 5), plot.title = element_blank(), axis.text = element_blank(),
          legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"), legend.text = element_text(size = 6),
          legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.83, 0.5), axis.title.y = element_text(size = 8))
  cluster_plots1[[spp]] = cluster_plots[[spp]] + annotation_map(shore, color = "black", fill = "white", size=0.1) +
    theme(text = element_text(size = 8), plot.title = element_blank())
  #mean_dens_plots[[spp]]$data$mean_exp <- log(mean_dens_plots[[spp]]$data$mean_exp) # switch back to log scale for visualization
  mean_dens_plots1[[spp]] = mean_dens_plots[[spp]] + annotation_map(shore, color = "black", fill = "white", size=0.1) +
    #scale_fill_viridis_c(trans="log",
    #                     breaks = c(min(mean_dens_plots[[spp]]$data$mean_exp), floor(max(mean_dens_plots[[spp]]$data$mean_exp)/200), floor(max(mean_dens_plots[[spp]]$data$mean_exp))),
    #                     labels = c(round(min(mean_dens_plots[[spp]]$data$mean_exp)), round(max(mean_dens_plots[[spp]]$data$mean_exp)/200, digits = -1), round(max(mean_dens_plots[[spp]]$data$mean_exp), digits = -2))) +
    theme(text = element_text(size = 8), plot.title = element_blank(),
          legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
          legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.83, 0.5))
  COG_plots1[[spp]] = COG_plots[[spp]] +
    scale_y_continuous(position = "right", breaks = seq(350, 500, by=50)) +
    coord_cartesian(ylim = c(355,539)) +
    theme(text = element_text(size = 8), plot.margin = unit(c(2.3,5.5,5.5,5.5), "pt"), legend.position = "none",
          plot.title = element_blank(), axis.text.x = element_blank()) +
    labs(tag = element_blank())
}

ggsave(filename = paste0("plots/Fig3_maps_COG_",knots,"_land.pdf"),
       plot = arrangeGrob(grobs = c(trend_plots1[spp_subset1],
                                    cluster_plots1[spp_subset1],
                                    mean_dens_plots1[spp_subset1],
                                    COG_plots1[spp_subset1]),
                          ncol = 4, as.table = FALSE, bottom = "Eastings (10s km)",
                          left = grid::textGrob("Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 6, height = 8.5, units = c("in"))

# For appendix figure verion of Fig 3, including all remaining species
ggsave(filename = paste0("plots/Fig3_maps_COG_",knots,"_land_SI.pdf"),
       plot = arrangeGrob(grobs = c(trend_plots1[spp_subset2],
                                    cluster_plots1[spp_subset2],
                                    mean_dens_plots1[spp_subset2],
                                    COG_plots1[spp_subset2]),
                          ncol = 4, as.table = FALSE, bottom = "Eastings (10s km)",
                          left = grid::textGrob("Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 6, height = 18, units = c("in"))

theme_set(theme_grey())
