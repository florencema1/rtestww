# ------------------------------------------
# Candidate number: 211392
# ------------------------------------------
# title: Deconvolution model outputs
# adapted from Wastewater-Based Estimation of the Effective Reproductive Number of SARS-CoV-2
# original author: Jana S Huisman et al.
# https://doi.org/10.1016/S2542-5196(21)00230-8
# repo: https://github.com/scwatson812/COVID19WastewaterModel
# ------------------------------------------

# Import EpiEstim for required functions
library(EpiEstim)

# Read in RNA level and case data
data <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Linked Data.csv")
case <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/cleaned cases.csv")
data$X <- NULL

# Read subsets into a list of dataframes
temp <- list.files(path="/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/subsets",
                  recursive = TRUE,
                  pattern = "\\.csv$",
                  full.names = TRUE)
myfiles <- lapply(temp, read.csv)

data$date <- lubridate::ymd(data$date)

d <- lapply(myfiles, function(df) {df["date"] <- lapply(df["date"], as.Date); df})

# Set empty dataframes for results
deconv_result <- list(data.frame())
result <- list(data.frame())

###########################################################
##### Deconvolve and Estimate WW Re #####

config_df = expand.grid("region" = c('south west', 'south east', 'london',
                                     'east of england', 'east midlands', 'west midlands',
                                     'north west', 'north east', 'yorkshire and the humber'),  
                        'incidence_var' = c('log_gc_per_day'),
                        'FirstGamma' = 'incubation',
                        'SecondGamma' = 'benefield' )


deconv_ww_data <- data.frame()
Re_ww <- data.frame()

for(row_i in 1:nrow(config_df)){
  new_deconv_data = deconvolveIncidence(data %>% filter(region == config_df[row_i, 'region']), 
                                        incidence_var = config_df[row_i, 'incidence_var'],
                                        getCountParams(as.character(config_df[row_i, 'FirstGamma'])), 
                                        getCountParams(as.character(config_df[row_i, 'SecondGamma'])),
                                        smooth_param = TRUE, n_boot = 50) #n_boot in paper: 1000
  
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

### Deconvolution and Re estimation for Case data ####
deconv_cases <- data.frame()
Re_cases <- data.frame()

config_df = expand.grid("region" = c('south west', 'south east', 'london',
                                     'east of england', 'east midlands', 'west midlands',
                                     'north west', 'north east', 'yorkshire and the humber'),  
                        'incidence_var' = c('newCasesBySpecimenDate'),
                        'FirstGamma' = 'incubation',
                        'SecondGamma' = 'benefield' )

for(inc_var in c('newCasesBySpecimenDate')){
  new_deconv_data = deconvolveIncidence(case, 
                                        incidence_var = inc_var,
                                        getCountParams('incubation'), 
                                        getCountParams(inc_var),
                                        smooth_param = TRUE, n_boot = 50) #n_boot = 1000 for paper
  new_Re = getReBootstrap(new_deconv_data) 
  
  deconv_cases <- bind_rows(deconv_cases, new_deconv_data)
  Re_cases = bind_rows(Re_cases, new_Re)
}

# Function
# for (j in c(1:9)) {
# for (incidence_var_i in c('log_gc_per_day')){
#   new_deconv_result[[j]] = deconvolveIncidence(d[[j]], 
#                                           incidence_var = incidence_var_i,
#                                           getCountParams('incubation'), 
#                                           getCountParams('benefield'),
#                                           smooth_param = TRUE, n_boot = 50) %>% #n_boot = 1000 in paper
#     mutate(data_type = incidence_var_i)
#   
#   new_result[[j]] = getReBootstrap(new_deconv_result[[j]])
#   new_result[[j]] = new_result[[j]] %>%
#     mutate(data_type = incidence_var_i)
#   new_result[[j]]['variable'] = incidence_var_i
#   
#   deconv_result[[j]] = bind_rows(deconv_result[[j]], new_deconv_result[[j]])
#   result[[j]] = bind_rows(result[[j]], new_result[[j]])
#   
# }
# }

# Plot bootstrapped result of Re
ggplot() +
  geom_ribbon(data=Re_ww, 
              aes(x=date, 
                  ymin = median_R_lowHPD,  
                  ymax = median_R_highHPD, 
                  color= region,
                  fill = data_type)) +
  geom_hline(yintercept = 1) +
  labs(title="Re Estimated by Deconvolution Method") +
  ylab("Re") +
  theme(legend.position = "none")

# Plot estimated number of cases
ggplot() +
  geom_point(data=deconv_ww_data,
             aes(x=date, y=value, color = region)) +
  labs(title="Incidence Estimated by Deconvolution Method") +
  ylab("Incidence") +
  theme(legend.position = "none")

## Deconvolution Plot ####
mean_deconv_data <- deconv_ww_data %>%
  group_by(date, region, country, source, data_type, incidence_var, 
           incubationParams, onsetToCountParams, GammaParams) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

deconv_ww_data <- data.frame()
Re_ww <- data.frame()

# deconv_plot <- ggplot() +
#   geom_errorbar(data = mean_deconv_data %>% filter(incidence_var == 'norm_n1'), 
#                 aes(x=date, ymin = value -sd,  ymax = value +sd, colour = GammaParams),
#                 show.legend = F, size = 1.2) +
#   labs(x = 'Date' , y='Estimated infection incidence') +
#   scale_x_date(date_breaks = "4 weeks", 
#                date_labels = '%b\n%d',
#                limits = c(as_date('2020-08-15'), as_date('2021-01-20'))
#   ) +
#   scale_colour_manual(values = viridis(4),
#                       breaks = c("incubation_zero","incubation_han",
#                                  "incubation_benefield", "incubation_death"),
#                       labels = c("Incubation only", "Incubation + Han", "Incubation + Benefield", "Incubation + Death"),
#                       name = 'Variable') +
#   theme(
#     strip.text.x= element_blank(),
#     axis.title.x =  element_blank()
#   )
# 
