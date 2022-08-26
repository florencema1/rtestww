# ----------------------------------------------------------------
# title: Exploratory data analysis of wastewater data of SARS-CoV-2 in from
#        45 sewage treatment works in England from July 2020 to February 2021
# data source: Morvan et al., UKHSA
# https://doi.org/10.1038/s41467-022-31753-y
# ----------------------------------------------------------------

library(tidyverse)
library(patchwork)

# Set themes for plots
theme_set(theme_minimal())

# Read in data
data <- read.csv("./data/linked_data.csv")
case <- read.csv("./data/cleaned_cases.csv")
r <- read.csv("./data/weighted_r.csv")
raw <- read.csv("./data/agg_data_with_region.csv")

# Remove index column 'X' if present
lst1 <- list(data=data, r=r, raw=raw)
list2env(lapply(lst1,`[`,-1), envir=.GlobalEnv)

# Convert date from string to date
data$date <- lubridate::ymd(data$date)
case$date <- lubridate::ymd(case$date)
r$date <- lubridate::ymd(r$date)
raw$date <- lubridate::ymd(raw$date)

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
raw <- arrange(mutate(raw,
                    region=factor(region,levels=neworder)),region)


# Weighted and interpolated VS actual RNA levels
g <- ggplot() + 
  geom_line(data, mapping = aes(x=date, y=weighted_avg_sars_cov2,
                                color='Interpolated')) +
  geom_line(data, mapping = aes(x=date, y=rna_7d_avg,
                                color='7-day Average')) +
  geom_point(raw, mapping = aes(x=date, y=sars_cov2_gc_l_mean), 
             colour="grey70", shape=1, size=0.1, alpha=0.5) +
  facet_wrap(~region,scales="free") +
  scale_color_manual(labels = c("7-day Average", "Interpolated"), 
                     values = c("tomato", "dodgerblue")) +
  labs(title="SARS-CoV-2 RNA Concentration - Weighted Average",
       x="Date (2020-2021)", 
       y="RNA Concentration (gc/l)",
       color="Legend") +
  coord_cartesian(ylim = c(0,5e5))
g + guides(alpha='none')

# Export plot
ggsave("plots/RNA over time.png", width = 8, height = 5, units = "in")

# Incidence over time
ggplot(case) + 
  geom_line(aes(x=date, y=newCasesBySpecimenDate, color='blue')) +
  geom_line(aes(x=date, y=cases_7d_avg, color='red')) + 
  facet_wrap(~region,scales="free") +
  labs(title = "Incidence of SARS-CoV-2", 
       x = "Date (2020-2021)", 
       y = "New Cases by Specimen Date", 
       color = "Legend") +
  scale_color_manual(labels = c("Actual case no.", "7-day average"), 
                     values = c("dodgerblue", "tomato")) +
  coord_cartesian(ylim = c(0,16000), xlim = c(as.Date("2020-09-04"),as.Date("2021-02-10")))

# Export plot 
ggsave("plots/Incidence over time.png", width = 8, height = 5, units = "in")

# Rt over time
g <- ggplot(r, aes(x=date, y=avg_r_avg)) + 
  geom_line() + 
  geom_ribbon(aes(x=date, 
                  ymin = avg_r_lower,
                  ymax = avg_r_upper, 
                  fill = NULL,
                  alpha = 0.3)) +
  geom_hline(yintercept=1, linetype="dashed") +
  facet_wrap(~region,scales="free") +
  labs(title="Effective Reproduction Number (Rt) - Weighted Average",
       x="Date (2020-2021)", y="Rt") +
  coord_cartesian(ylim = c(0.6,1.67), xlim = c(as.Date("2020-09-04"),as.Date("2021-02-10"))) +
  theme(legend.position = "none") 
g + guides(alpha="none")

# Export plot 
ggsave("plots/Rt over time.png", width = 8, height = 5, units = "in")


# # Trial to combine RNA and incidence plots (2 axis scales)
# max_rna <- max(data$weighted_avg_sars_cov2)
# max_rolling <- max(case$rolling)
# 
# d2 <- gather(d1, 'var', 'val', stones:revenue) %>% 
#   mutate(val = if_else(var == 'revenue', as.double(val), val / (max_rolling / max_rna)))
# 
# ggplot() + 
#   geom_line(data, mapping = aes(x=date, y=gene_copies_per_day)) +
#   geom_line(case, mapping = aes(x=date, y=rolling)) +
#   geom_point(raw, mapping = aes(x=date, y=sars_cov2_gc_l_mean), 
#              colour="red", size=0.1, alpha=0.3) +
#   facet_wrap(~region,scales="free") +
#   labs(title="Weighted and Interpolated Level of SARS-CoV-2") +
#   xlab("Date (2020-2021)") +
#   ylab("Weighted Average RNA Level (gc/l)") +
#   scale_y_continuous(limits=c(0,4e5), sec.axis = sec_axis(trans = ~ . * (max_stones / max_revenue),
#                                                           name = 'number of stones'),
#                      labels = dollar)
