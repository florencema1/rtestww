# ----------------------------------------------------------------
# Title: Estimations from deconvolution model
# Functions from: Wastewater-Based Estimation of the Effective 
# Reproductive Number of SARS-CoV-2 by Jana S Huisman et al.
# https://doi.org/10.1016/S2542-5196(21)00230-8
# Repo: https://github.com/scwatson812/COVID19WastewaterModel
# ----------------------------------------------------------------

# Import EpiEstim for required functions
library(tidyverse)
library(EpiEstim)
library(patchwork)

# Set themes for plots
theme_set(theme_minimal())

# Read in data
data <- read.csv("./data/linked_data.csv")
case <- read.csv("./data/cleaned_cases.csv")

# Convert date from string to date
data$date <- lubridate::ymd(data$date)
case$date <- lubridate::ymd(case$date)

# Re-order groups for plots
neworder <- c('North West', 'Yorkshire And The Humber', 'North East',
              'West Midlands', 'East Midlands', 'East Of England',
              'South West','London' , 'South East')
data <- arrange(mutate(data,
                       region=factor(region,levels=neworder)),region)
case <- arrange(mutate(case,
                       region=factor(region,levels=neworder)),region)
r <- arrange(mutate(r,
                    region=factor(region,levels=neworder)),region)

###########################################################
## Normalisation of WW data ####

# Min gene copies per day observed for each region
norm_min <- data %>%
  group_by(region) %>%
  summarise(norm_factor = min(gc_per_day))

write.csv(norm_min, 'normalisation_factor.csv')

# Normalisation
ww_data <- data %>%
  merge(norm_min, by='region') %>%
  mutate(norm_gc = gc_per_day/norm_factor)

# # (constant normalisation factor used in the original study)
# ww_data <- data %>%
#   mutate(norm_gc = gc_per_day/1e12)

# Check normalised levels across different regions
# Colors based on median value
wd <- ww_data %>% mutate(region=region) %>%
  group_by(region) %>%
  mutate(med_ngc = median(norm_gc))

# Quantiles required for boxplot
custom_quantile <- function(x) {
  out <- c(quantile(x, 0), quantile(x, 0.25), median(x), quantile(x,0.75), quantile(x, 1))
  names(out) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(out)
}

# Re-order groups
wd <- arrange(mutate(wd,region=factor(region,levels=neworder)),region)

ggplot(wd, aes(x=factor(region), y=norm_gc, fill=med_ngc)) + 
  stat_summary(fun.data=custom_quantile, geom='boxplot') +
  labs(title='Normalised RNA Levels',
       x='Region',
       y='RNA Level (gene copies per day)',
       fill='Median') +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 20, hjust=1)) 
  
ggsave("plots/Normalisation Deconvolution.png", width = 8, height = 4.5, units = "in")

############################################################
### Deconvolution and Re estimation for Wastewater data ####

config_df = expand.grid("region" = c('North West', 'Yorkshire And The Humber', 'North East',
                                     'West Midlands', 'East Midlands', 'East Of England',
                                     'South West','London' , 'South East'),  
                        'incidence_var' = 'norm_gc',
                        'FirstGamma' = 'incubation',
                        'SecondGamma' = 'benefield' )


deconv_ww_data <- data.frame()
Re_ww <- data.frame()

for(row_i in 1:nrow(config_df)){
  new_deconv_data = deconvolveIncidence(ww_data %>% filter(region == config_df[row_i, 'region']), 
                                        incidence_var = config_df[row_i, 'incidence_var'],
                                        getCountParams(as.character(config_df[row_i, 'FirstGamma'])), 
                                        getCountParams(as.character(config_df[row_i, 'SecondGamma'])),
                                        smooth_param = TRUE, n_boot = 50) 
  
  new_deconv_data <- new_deconv_data %>%
    mutate(incidence_var = config_df[row_i, 'incidence_var'])
  
  ##### Get Re #####
  new_Re_ww = getReBootstrap(new_deconv_data)
  new_Re_ww <- new_Re_ww %>%
    mutate(variable = config_df[row_i, 'incidence_var'],
           region = config_df[row_i, 'region'])
  
  deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
  Re_ww = bind_rows(Re_ww, new_Re_ww)
}

# Plot bootstrapped result of Re
g <- ggplot() + 
  geom_ribbon(data=Re_ww, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4)) +
  geom_line(data=Re_ww, aes(x = date, y = median_R_mean), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_hline(yintercept=1, linetype='dashed') +
  facet_wrap(~region,scales="free") +
  labs(title="Rt Estimated from Wastewater Data",
       x="Date (2020-2021)",
       y="Rt") +
  coord_cartesian(ylim = c(0, 3)) 
g + scale_alpha(guide = 'none')

# Export plot and results
ggsave("plots/Deconvolution R WW.png", width = 8, height = 4.5, units = "in")
write_csv(Re_ww, "Deconvolution R from WW.csv")

# Prepare data for incidence plot
plot_data <- deconv_ww_data %>%
  select(-local_infection, -country, -source) %>%
  group_by(date, region, data_type) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

# Re-order groups
plot_data <- arrange(mutate(plot_data,
                            region=factor(region,levels=neworder)),region)

# Plot estimated incidence
ax <- ggplot() +
  geom_ribbon(data=plot_data,
             aes(x=date,
                 y=value,
                 ymin=value-sd,
                 ymax=value+sd),
                 alpha=0.5) +
  geom_line(data=plot_data,
            aes(x = date, y = value), 
            alpha = 1, size = 0.5, show.legend=F) +
  labs(title="Incidence Estimated by Deconvolution Method") +
  ylab("Incidence (cases/day)") +
  xlab("Date (2020-2021)") +
  facet_wrap(~region,scales="free") 

ax

ggsave("plots/Deconvolution Incidence.png", width = 8, height = 4.5, units = "in")

# Incidence over time
bx <- ggplot(case) + 
  geom_line(aes(x=date, y=cases_7d_avg)) + 
  facet_wrap(~region,scales="free") +
  labs(title = "New Cases by Specimen Date", 
       x = "Date (2020-2021)", y = "Incidence") +
  coord_cartesian(xlim=c(as.Date("2020-09-04"),as.Date("2021-02-10"))) 

ax/bx

# Export plot and results
ggsave("plots/Deconvolution Incidence (Combined).png", width = 8, height = 10, units = "in")
write_csv(Re_ww, "Deconvolution Incidence from WW.csv")

######################################################
### Deconvolution and Re estimation for Case data ####
config_case = expand.grid("region" = c('North West', 'Yorkshire And The Humber', 'North East',
                                       'West Midlands', 'East Midlands', 'East Of England',
                                       'South West','London' , 'South East'),  
                        'incidence_var' = 'newCasesBySpecimenDate',
                        'FirstGamma' = 'incubation',
                        'SecondGamma' = 'zero' )


deconv_cases <- data.frame()
Re_cases <- data.frame()
for(row_i in 1:nrow(config_case)){
  new_deconv_data = deconvolveIncidence(case %>% filter(region == config_case[row_i, 'region']), 
                                        incidence_var = config_case[row_i, 'incidence_var'],
                                        getCountParams(as.character(config_case[row_i, 'FirstGamma'])), 
                                        getCountParams(as.character(config_case[row_i, 'SecondGamma'])),
                                        smooth_param = TRUE, n_boot = 50) 
  
  new_deconv_data <- new_deconv_data %>%
    mutate(incidence_var = config_df[row_i, 'incidence_var'])
  
  ##### Get Re #####
  new_Re = getReBootstrap(new_deconv_data)
  new_Re <- new_Re %>%
    mutate(variable = config_case[row_i, 'incidence_var'],
           region = config_case[row_i, 'region'])
  
  deconv_cases <- bind_rows(deconv_cases, new_deconv_data)
  Re_cases = bind_rows(Re_cases, new_Re)
}

# Plot bootstrapped result of Re
g <- ggplot() + 
  geom_ribbon(data=Re_cases, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4)) +
  geom_line(data=Re_cases, aes(x = date, y = median_R_mean), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_hline(yintercept=1, linetype='dashed') +
  facet_wrap(~region,scales="free") +
  labs(title="Rt Estimated from Case Data") +
  xlab("Date (2020-2021)") +
  ylab("Rt") +
  coord_cartesian(ylim = c(0.5, 2)) 
g + scale_alpha(guide = 'none')

# Export plot and results
ggsave("plots/Deconvolution R Cases.png", width = 8, height = 4.5, units = "in")
write_csv(Re_ww, "Deconvolution R from Cases.csv")

###########################################################
### RMSE for deconvolution model (wastewater vs case data) 
### and deconvolution vs Reported Rt (wastewater vs report)

rmse_all <- data.frame()

for(row_i in 1:nrow(config_case)){
  new_se <- Re_ww %>% filter(region == config_case[row_i, 'region']) %>%
  left_join(Re_cases, by = c('date', 'region', 'country', 'source'), 
            suffix = c('.i', '.j')) %>%
  left_join(data, by= c('date', 'region')) %>%
  mutate(se_est = (median_R_mean.i - median_R_mean.j)^2) %>%
  mutate(se_report = (median_R_mean.i - avg_r_avg)^2)
  
  new_rmse <- cbind(rmse_est = sqrt(sum(new_se$se_est, na.rm = T)/length(new_se$date)),
                    rmse_report = sqrt(sum(new_se$se_report, na.rm = T)/length(new_se$date)))
  
  rmse_all <- rbind(rmse_all, new_rmse)
}

rmse_all <- cbind(region=config_case$region, rmse_all)
write.csv(rmse_all, 'RMSE.csv')

###########################################################
### Combined plots ###

# Rt
g <- ggplot() + 
  geom_ribbon(data=Re_cases, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4,
                  color="Deconvolution\n(Case data)",
                  fill="Deconvolution\n(Case data)")) +
  geom_line(data=Re_cases, aes(x = date, y = median_R_mean, 
                               color="Deconvolution\n(Case data)"), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_ribbon(data=Re_ww, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4,
                  color="Deconvolution\n(Wastewater data)",
                  fill="Deconvolution\n(Wastewater data)")) +
  geom_line(data=Re_ww, aes(x = date, y = median_R_mean, 
                            color="Deconvolution\n(Wastewater data)"), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_ribbon(data=data, aes(x=date, 
                  ymin = avg_r_lower,
                  ymax = avg_r_upper, 
                  color="Reported Rt",
                  fill = "Reported Rt",
                  alpha = 0.4)) +
  geom_line(data=r, aes(x=date, y=avg_r_avg, color="Reported Rt")) +
  geom_hline(yintercept=1, linetype='dashed') +
  facet_wrap(~region,scales="free") +
  labs(title="Rt Estimates from Deconvolution Model and Reported Rt", 
       x="Date (2020-2021)", y="Rt", fill='Rt') +
  coord_cartesian(ylim = c(0,3), xlim = c(as.Date("2020-09-04"),as.Date("2021-02-05"))) 
g + guides(alpha="none", colour='none')

ggsave("plots/Deconvolution R (All).png", width = 8, height = 4.5, units = "in")

g <- ggplot() + 
  geom_ribbon(data=Re_cases, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4,
                  color="Deconvolution\n(Case data)",
                  fill="Deconvolution\n(Case data)")) +
  geom_line(data=Re_cases, aes(x = date, y = median_R_mean, 
                               color="Deconvolution\n(Case data)"), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_ribbon(data=data, aes(x=date, 
                          ymin = avg_r_lower,
                          ymax = avg_r_upper, 
                          color="Reported Rt",
                          fill = "Reported Rt",
                          alpha = 0.4)) +
  geom_line(data=r, aes(x=date, y=avg_r_avg, color="Reported Rt")) +
  geom_hline(yintercept=1, linetype='dashed') +
  facet_wrap(~region,scales="free") +
  labs(title="Rt Estimates from Deconvolution Model (Case Data) and Reported Rt", 
       x="Date (2020-2021)", y="Rt", fill='Rt') +
  coord_cartesian(ylim = c(0.6,1.67), xlim = c(as.Date("2020-09-04"),as.Date("2021-02-05"))) 
g + guides(alpha="none", colour='none')

ggsave("plots/Deconvolution R (case vs report).png", width = 8, height = 4.5, units = "in")

g <- ggplot() + 
  geom_ribbon(data=Re_cases, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4,
                  color="Deconvolution\n(Case data)",
                  fill="Deconvolution\n(Case data)")) +
  geom_line(data=Re_cases, aes(x = date, y = median_R_mean, 
                               color="Deconvolution\n(Case data)"), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_ribbon(data=Re_ww, 
              aes(x=date, 
                  ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, 
                  alpha=0.4,
                  color="Deconvolution\n(Wastewater data)",
                  fill="Deconvolution\n(Wastewater data)")) +
  geom_line(data=Re_ww, aes(x = date, y = median_R_mean, 
                            color="Deconvolution\n(Wastewater data)"), 
            alpha = 1, size = 0.5, show.legend=F) +
  geom_hline(yintercept=1, linetype='dashed') +
  facet_wrap(~region,scales="free") +
  labs(title="Rt Estimates from Deconvolution Model", 
       x="Date (2020-2021)", y="Rt", fill='Rt') +
  coord_cartesian(ylim = c(0,3), xlim = c(as.Date("2020-09-04"),as.Date("2021-02-05")))
g + guides(alpha="none", colour="none")

ggsave("plots/Deconvolution R (ww vs case).png", width = 8, height = 4.5, units = "in")

