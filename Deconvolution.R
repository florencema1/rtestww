### Deconvolution of wastewater data by JS Huisman et al. ###
###### Deconvolve #####

# Read in RNA level
data <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Linked Data_South West.csv")
data$X <- NULL

# Read subsets into a list of dataframes
temp <- list.files(path="/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/subsets",
                  recursive = TRUE,
                  pattern = "\\.csv$",
                  full.names = TRUE)
myfiles <- lapply(temp, read.csv)

data$date <- lubridate::ymd(data$date)

# Set empty dataframes for results
deconv_result <- data.frame()
result <- data.frame()

# Function
for (incidence_var_i in c('log_gene')){
  new_deconv_result = deconvolveIncidence(data, 
                                          incidence_var = incidence_var_i,
                                          getCountParams('incubation'), 
                                          getCountParams('benefield'),
                                          smooth_param = TRUE, n_boot = 50) %>% #n_boot = 1000 in paper
    mutate(data_type = incidence_var_i)
  
  new_result = getReBootstrap(new_deconv_result)
  new_result = new_result %>%
    mutate(data_type = incidence_var_i)
  new_result['variable'] = incidence_var_i
  
  deconv_result = bind_rows(deconv_result, new_deconv_result)
  result = bind_rows(result, new_result)
  
}

# Plot bootstrapped result of Re
ggplot() +
  geom_ribbon(data=result, 
              aes(x=date, ymin = median_R_lowHPD,  ymax = median_R_highHPD, fill = data_type)) +
  labs(title="Re Estimated by Deconvolution Method (SW)") +
  ylab("Re") +
  theme(legend.position = "none")

# Plot estimated number of cases
ggplot() +
  geom_point(data=deconv_result,
             aes(x=date, y=value, color = data_type)) +
  labs(title="Incidence Estimated by Deconvolution Method (SW)") +
  ylab("Incidence") +
  theme(legend.position = "none")

