library(sdmTMB)
library(ggplot2)
library(dplyr)
dir.create("figs", showWarnings = FALSE)

set.seed(122)
d <- expand.grid(X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100), year = 1:5)

.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}

sigma_O <- 0.3
sigma_Z <- 0.3
sigma_E <- 0.15
kappa <- 1.7
rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1/kappa)
d$omega_s <- .rf_sim(model = rf_omega, d$X, d$Y)

# ggplot(dplyr::filter(d, year == 1), aes(X, Y, fill = omega_s)) +
#   geom_raster() +
#   scale_fill_gradient2()

rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1/kappa)

d <- d %>% group_by(year) %>%
  group_split() %>%
  purrr::map_dfr(
    ~tibble(X = .x$X, Y = .x$Y, omega_s = .x$omega_s, year = .x$year,
      eps_st = .rf_sim(model = rf_epsilon, .x$X, .x$Y))) %>%
  ungroup()

# ggplot(d, aes(X, Y, fill = eps_st)) +
#   geom_raster() +
#   scale_fill_gradient2() +
#   facet_grid(cols = vars(year))
#
# ggplot(d, aes(X, Y, fill = eps_st + omega_s)) +
#   geom_raster() +
#   scale_fill_gradient2() +
#   facet_grid(cols = vars(year))

rf_zeta <- RandomFields::RMmatern(nu = 1, var = sigma_Z^2, scale = 1/kappa)
d$zeta_s <- .rf_sim(model = rf_zeta, d$X, d$Y)

d$year_cent <- d$year - mean(d$year)
d$`With spatially\nvarying trend` <- d$omega_s + d$eps_st + d$zeta_s * d$year_cent
d$`Without spatially\nvarying trend` <- (d$omega_s + d$eps_st) * 1.5 # approx. equalize variance

g <- reshape2::melt(d, id.vars = c("X", "Y", "year"),
  measure.vars = c("With spatially\nvarying trend", "Without spatially\nvarying trend")) %>%
  ggplot(aes(X, Y, fill = value)) +
  geom_raster() +
  # scale_fill_gradient2(high = scales::muted("red"), mid = "grey94",
  #   low = scales::muted("blue")) +
  scale_fill_viridis_c(option = "A") +
  facet_grid(cols = vars(year), rows = vars(variable)) +
  ggsidekick::theme_sleek() +
  theme(panel.border = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
    panel.spacing.x = unit(-0.15, "lines"), panel.spacing.y = unit(-0.25, "lines"),
    axis.ticks = element_blank()) +
  guides(fill = FALSE)

# print(g)
ggsave("figs/fig1.pdf", width = 5, height = 2.5)

