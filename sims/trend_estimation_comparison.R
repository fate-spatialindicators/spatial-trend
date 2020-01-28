library(dplyr)
library(ggplot2)
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
set.seed(123)

# create time series with varing slopes
n = 1000 # total datasets
t = 15 # time series length
# generate data where there's at least one obs/year, but a couple obs in some years
df = expand.grid("sim"=1:n, "time"=1:t)
df = rbind(df, expand.grid("sim"=1:n, "time"=seq(1,t,by=2)))

# generate low/medium/high estimates of trends,
# and low/med/high variances of random effects
pars = data.frame("sim"=1:n,
  "b0"=rnorm(n),
  "b1"=sample(c(0.1,0.2,0.4),size=n,replace=T),
  "re"=sample(c(0.01,0.05,0.1),size=n,replace=T))
# join in design and parameters
df = left_join(df,pars)

# calculate predictions: E[y] = B0 + B1*t + re(t)
df$pred = df$b0 + df$b1*df$time + rnorm(nrow(df),0,sd=df$re)
# add observation error -- could be confounded with re(t) if one obs/year
df$y = rnorm(nrow(df), df$pred, 0.1)
pars$post_hoc_trend=NA
pars$glmm_trend=NA

for(i in 1:n) {
  # fit model without trend and do regression post-hoc
  mod = glmmTMB::glmmTMB(y ~ (1|time), data = dplyr::filter(df,sim==i))
  pars$post_hoc_trend[i] = coef( lm(coef(mod)$cond$time$`(Intercept)` ~ seq(1,t)))[2]

  # include trend estimate directly
  mod = glmmTMB::glmmTMB(y ~ time + (1|time), data = dplyr::filter(df,sim==i))
  pars$glmm_trend[i] = coef(mod)$cond$time$time[1]
  print(i)
}

saveRDS(pars,file="plots/trend_comparison.rds")


g1 = ggplot(pars, aes(as.factor(b1),post_hoc_trend-b1)) +
  geom_hline(aes(yintercept=0),col="grey50") +
  geom_boxplot(fill="dark blue",alpha=0.3,col="dark blue") +
  xlab("B[1]") +
  ylab(expression(paste("Trend bias (",hat(theta)," - ",theta,")"))) +
  ggtitle("Post-hoc trend estimation") +
  theme_sleek()

g2 = ggplot(pars, aes(as.factor(re),post_hoc_trend-b1)) +
  geom_hline(aes(yintercept=0),col="grey50") +
  geom_boxplot(fill="dark blue",alpha=0.3,col="dark blue") +
  xlab(expression(paste("Random effect ",sigma))) +
  ylab(expression(paste("Trend bias (",hat(theta)," - ",theta,")"))) +
  ggtitle("Post-hoc trend estimation")

g3 = ggplot(pars, aes(as.factor(b1),glmm_trend-b1)) +
  geom_hline(aes(yintercept=0),col="grey50") +
  geom_boxplot(fill="dark blue",alpha=0.3,col="dark blue") +
  xlab("B[1]") +
  ylab(expression(paste("Trend bias (",hat(theta)," - ",theta,")"))) +
  ggtitle("Trend estimation within GLMM")

g4 = ggplot(pars, aes(as.factor(re),glmm_trend-b1)) +
  geom_hline(aes(yintercept=0),col="grey50") +
  geom_boxplot(fill="dark blue",alpha=0.3,col="dark blue") +
  geom_hline(aes(yintercept=0)) +
  xlab(expression(paste("Random effect ",sigma))) +
  ylab(expression(paste("Trend bias (",hat(theta)," - ",theta,")"))) +
  ggtitle("Trend estimation within GLMM")

pdf("plots/Figure S1 trend_estimation_comparison.pdf")
gridExtra::grid.arrange(g1, g2, g3, g4, nrow=2)
dev.off()
