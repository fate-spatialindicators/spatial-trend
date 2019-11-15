library(sdmTMB)
library(ggplot2)
library(viridis)
library(dplyr)
library(fpc)
library(gridExtra)

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

# define species groups, in alphabetic order of species
species_groups = c("flatfish","chond","rockfish","rockfish","rockfish",
                  "flatfish", "flatfish", "roundfish", "chond",
                  "flatfish", "rockfish", "flatfish", "flatfish", "roundfish",
                  "rockfish", "chond", "rockfish", "chond", "rockfish")

anisotropy_plots = list()
qq_plots = list()
residuals_plots = list()
prediction_plots = list()
prediction_plots_ex = list()
spatiotemporal_plots = list()
trend_plots = list()
intercept_plots = list()
spatial_plots = list()
cluster_plots = list()
COG_plots = list()
all_clust = NULL
st_CVs = NULL

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
     p = predict(d, newdata=wc_grid_all)

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ## model checking
     # check anisotropy
     #anisotropy_plots[[spp]] = plot_anisotropy(d) +
     # ggtitle(paste0(species[spp],"_aniso"))
     ## residuals
     #data = select(d$data, X, Y, cpue_kg_km2, year)
     #data$residuals = residuals(d)
     # qq plots
     #qq_plots[[spp]] = qqnorm(data$residuals)
     # spatial residuals
     #residuals_plots[[spp]] = plot_map_point(data, "residuals") + facet_wrap(~year) + geom_point(size=0.05, alpha=0.1) +
      # coord_fixed() + scale_color_gradient2() + ggtitle(species[spp])
     # check convergence
     #sd = as.data.frame(summary(TMB::sdreport(d$tmb_obj)))
     #print(species[spp])
     #print(d$sd_report)
     # check whether AR1 assumption is supported in models where fields are not IID, printing estimate and 95%CI for AR1 param
     #print("AR1 estimate")
     #print(sd$Estimate[row.names(sd) == "ar1_phi"])
     #print("AR1 parameter 95% CI")
     #print(sd$Estimate[row.names(sd) == "ar1_phi"] +
     # c(-2, 2) * sd$`Std. Error`[row.names(sd) == "ar1_phi"])

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     #maps
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
     # make plot of spatial intercept
     intercept_plots[[spp]] = plot_map_raster(dplyr::filter(p,year==2003), "omega_s") +
       scale_fill_gradient2() +
       theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
             axis.title = element_blank(),
             axis.text = element_blank(),
             legend.key.width = unit(0.2,"cm"),
             legend.title = element_blank(),
             legend.position = "right",
             legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
       labs(title = "Intercept")

     # make plot of predictions from full model for example year
     prediction_plots_ex[[spp]] = plot_map_raster(dplyr::filter(p,year==ex_year), "est") +
       scale_fill_viridis_c() +
       theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
             axis.title = element_blank(),
             axis.text = element_blank(),
             legend.key.width = unit(0.2,"cm"),
             legend.title = element_blank(),
             legend.position = "right",
             legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
       labs(title = "Density")

     ## plot absolute rate of change in space
     # extract all coefficients
     #sr_se <- summary(d$sd_report)[,"Std. Error"]
     #b_j <- unname(d$model$par[grep("b_j", names(d$model$par))])
     #b_j_se <- unname(sr_se[grep("b_j", names(sr_se))])
     #mm <- cbind(b_j, b_j_se)
     #colnames(mm) <- c("coef.est", "coef.se")
     #row.names(mm) <- colnames(model.matrix(d$formula, d$data))
     # get mean coefficient of year effects
     #coef_fixed <- mean(unname(d$model$par[grep("b_j", names(d$model$par))])[c(1,4:18)])
     # get absolute spatial change by combining mean of fixed effects with spatial trend, intercept
     # (TO DO: choose numeric time t to multiply by slope field, and factor to multiply by mean year coefficient)
     # most recent year would be 7 or 7.5 depending on how t_i handled. use 0 for log depth, 1 for log depth squared
     #p$abs <-  exp(dplyr::filter(p,year==2003)$omega_s + 7*dplyr::filter(p,year==2003)$zeta_s + z*coef_fixed)

     # show example of "static" spatial components to visually weight trend by adding the intercept random field to the fixed effects,
     # giving the predictions without the slope field and spatiotemporal fields
     p$omega_fixed <- p$omega_s + p$est_non_rf

     spatial_plots[[spp]] = plot_map_raster(dplyr::filter(p,year==ex_year), "omega_fixed") +
       scale_fill_viridis_c() +
       theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5),
             axis.title = element_blank(),
             axis.text = element_blank(),
             legend.key.width = unit(0.2,"cm"),
             legend.title = element_blank(),
             legend.position = "right",
             legend.margin=margin(c(0,0,0,-1.2), unit='cm')) +
       labs(title = "Trend weight")

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
     # Adding legend manually because in the code below, legend.position to put legend inside plot causes coastwide COG line to disapear!
     #scale_color_gradient(low = "#feb24c", high = "#a50f15", "Biomass", guide = "colourbar", breaks = c(0,1), limits = c(0, 1)) +
     #theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
     #     axis.text.x = element_blank(), plot.margin = unit(c(0,0,1,3), "pt"),
     #    legend.key.height = unit(0.5,"line"), legend.key.width = unit(0.8, "line"),
     #   legend.margin = margin(c(-1.9,1,0,-9.4), unit='line'), legend.direction = "horizontal",legend.justification = c(0,1),
     #  legend.title = element_text(size = 12, face = "bold")) +
     # guides(color = guide_colorbar(title.position="top", title.hjust = 0.5))

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # save spatial and temporal CVs for each species
     st_CVs = rbind(st_CVs,
                    cbind(p %>% filter(year==2003) %>% summarize(cv_s = sd(exp(omega_s))/mean(exp(omega_s))),
                    p %>% group_by(year) %>% summarize(mean_est = mean(exp(est))) %>% ungroup() %>%
                      summarize(cv_t = sd(mean_est)/mean(mean_est)),
                    species = species[spp]))
}

# save results
save.image(file = paste0("results/plotdata_lim_yearfixed_",knots,".Rdata"))
#save(vg_stats, ref_ranges, all_clust, st_CVs, st_CVs_cluster, file = "results/plotdata.Rdata")

# COG timeseries plots
ggsave(filename = paste0("plots/FigS2_COG_color",knots,".pdf"),
       plot = arrangeGrob(grobs = COG_plots, ncol = 5, bottom = "year",
                          left = grid::textGrob("COG Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 12, height = 8, units = c("in"))

# plot comparison of spatial and temporal CVs by species, species group
st_CVs$species_group = species_groups
ggsave(filename = paste0("plots/SpatialTemporalCV",knots,".pdf"),
       plot = ggplot(st_CVs, aes(x=cv_s, y=cv_t, color = species_group, label = species)) +
              geom_point() +
              labs(x = "Spatial CV", y = "Temporal CV"), # + geom_text(aes(label = species), size = 2, hjust=0.5, vjust=-0.6),
       width = 6, height = 4, units = c("in"))

# plot of predictions from full model (all fixed + random effects)
ggsave(filename = paste0("plots/predicted_density_maps",knots,".pdf"),
       plot = marrangeGrob(prediction_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))

# plot only spatiotemporal random effects
ggsave(filename = paste0("plots/st_maps_",knots,".pdf"),
       plot = marrangeGrob(spatiotemporal_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))

# plot COG along with trend, intercept predicted to regular grid, along with clusters, for selected species
spp_subset1 = c(1, 3, 7, 12, 14, 15)
trend_plots1 = list()
spatial_plots1 = list()
cluster_plots1 = list()
intercept_plots1 = list()
prediction_plots_ex1 = list()
COG_plots1 = list()

for(spp in 1:length(species)){
  if(spp == spp_subset1[1]){
    trend_plots1[[spp]] = trend_plots[[spp]] +
      theme(text = element_text(size = 5), plot.title = element_text(size=9, margin = margin(t = 5, b = -20), hjust = 0.55),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"), legend.text = element_text(size = 6),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5), axis.title.y = element_text(size = 8))
    cluster_plots1[[spp]] = cluster_plots[[spp]] +
      theme(text = element_text(size = 8), plot.title = element_text(size=9, margin = margin(t = 5, b = -20), hjust = 0.62))
    spatial_plots1[[spp]] = spatial_plots[[spp]]+
      theme(text = element_text(size = 8), plot.title = element_text(size=9, margin = margin(t = 5, b = -20), hjust = 0.55),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5))
    prediction_plots_ex1[[spp]] = prediction_plots_ex[[spp]]+
      theme(text = element_text(size = 8), plot.title = element_text(size=9, margin = margin(t = 5, b = -20), hjust = 0.55),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5))
    intercept_plots1[[spp]] = intercept_plots[[spp]] +
      theme(text = element_text(size = 8), plot.title = element_text(size=9, margin = margin(t = 5, b = -20), hjust = 0.55),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5))
    COG_plots1[[spp]] = COG_plots[[spp]] +
      scale_y_continuous(position = "right", breaks = seq(350, 500, by=50)) +
      coord_cartesian(ylim = c(355,539)) +
      theme(text = element_text(size = 8), plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"), legend.position = "none",
            plot.title = element_text(margin = margin(t = 5, b = -20), hjust = 0.55)) +
      labs(title = "COG", tag = element_blank())
  }
  else{
    trend_plots1[[spp]] = trend_plots[[spp]] +
      theme(text = element_text(size = 5), plot.title = element_blank(),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"), legend.text = element_text(size = 6),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5), axis.title.y = element_text(size = 8))
    cluster_plots1[[spp]] = cluster_plots[[spp]] +
      theme(text = element_text(size = 8), plot.title = element_blank())
    spatial_plots1[[spp]] = spatial_plots[[spp]] +
      theme(text = element_text(size = 8), plot.title = element_blank(),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5))
    prediction_plots_ex1[[spp]] = prediction_plots_ex[[spp]] +
      theme(text = element_text(size = 8), plot.title = element_blank(),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5))
    intercept_plots1[[spp]] = intercept_plots[[spp]] +
      theme(text = element_text(size = 8), plot.title = element_blank(),
            legend.key.height = unit(0.55, "line"), legend.key.width = unit(0.2, "line"),
            legend.margin=margin(c(0,0,0,0), unit='cm'), legend.position = c(0.93, 0.5))
    COG_plots1[[spp]] = COG_plots[[spp]] +
      scale_y_continuous(position = "right", breaks = seq(350, 500, by=50)) +
      coord_cartesian(ylim = c(355,539)) +
      theme(text = element_text(size = 8), plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"), legend.position = "none",
            plot.title = element_blank()) +
      labs(tag = element_blank())
  }
  # was trying to add this to above loop to make x axis labels on bottom panels,
  # but there is an issue with the previous theme settings being empty as the original plot did not have axes
  #if(spp == 19){
  #  trend_plots1[[spp]] = trend_plots1[[spp]] + theme(axis.text.x = element_text())
  #  trend_cluster_plots1[[spp]] = trend_cluster_plots1[[spp]] + theme(axis.text.x = element_text())
  #  intercept_plots1[[spp]] = intercept_plots1[[spp]] + theme(axis.text.x = element_text())
  #  intercept_cluster_plots1[[spp]] = intercept_cluster_plots1[[spp]] + theme(axis.text.x = element_text())
  #  COG_plots1[[spp]] = COG_plots1[[spp]] +
  #    theme(axis.text.x = element_text(), axis.title = element_text(size = 8),
  #          legend.key.height = unit(0.2, "line"), legend.key.width = unit(0.8, "line"),
  #          legend.direction = "horizontal", legend.title = element_text(size = 8, face = "bold")) +
  #    guides(color = guide_colorbar(title.position="bottom", title.hjust = 0.5))
  #}
}
ggsave(filename = paste0("plots/Fig3_maps_COG_",knots,".pdf"),
       plot = arrangeGrob(grobs = c(trend_plots1[spp_subset1],
                                    cluster_plots1[spp_subset1],
                                    spatial_plots1[spp_subset1],
                                    COG_plots1[spp_subset1]),
                          ncol = 4, as.table = FALSE, bottom = "Eastings (10s km)",
                          left = grid::textGrob("Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 6.5, height = 8.5, units = c("in"))

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
ggsave(filename = paste0("plots/Fig4_trendcluster_stripplot_",knots,".pdf"),
       width = 10, height = 6, units = c("in"))
