# ------------------------------------------
# title: Transforming RNA level into case number
# ------------------------------------------
library(tidyverse)
library(linelist)
library(EpiEstim)

# Set ggplot theme
theme_set(theme_bw())

# Read in RNA level
data <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Linked Data_South West.csv")
data$X <- NULL

# Transform into date format
data$date <- lubridate::ymd(data$date)
max(data$date)-min(data$date) 

# Read in estimated case numbers from COVID19WasteWater SEIR model
case <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/model_estimates.csv")
model_factor <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/incidence_variables.csv")

# Insert index number as 'day number' for joining tables later
df <- tibble::rowid_to_column(data, "day")
case <- tibble::rowid_to_column(case, "day")

# Estimate case number from RNA level
bx <- bs(combined$rna, df=10)
est_med_med <- cbind(1,bx) %*% table[[1]]$coef3 #coef 3 corresponds to median
est_med_low <- cbind(1,bx) %*% table[[1]]$coef1 #coef 1, 5 refer to 95% CI; coef 2, 4 refer to 75% CI
est_med_high <- cbind(1,bx) %*% table[[1]]$coef5 

# Create dataframe for plotting
sw_est <- as.data.frame(cbind(day = combined$day, 
                              rna  = combined$rna)) %>%
  cbind(activecase_m=est_med_med) %>%
  cbind(activecase_l=est_med_low) %>%
  cbind(activecase_h=est_med_high) 

# Export est. incidence 
write.csv(sw_est,"sw_est.csv", row.names=FALSE)

# Pivot dataframe for plotting multiple lines on the same graph
plotdf <- sw_est%>%
          pivot_longer('rna':'activecase_m')

# # Plot but both lines have ribbon
# ggplot(plotdf, aes(x=week, y=value, group=name, color=name)) + 
#   geom_line() +
#   labs(title='RNA level (gene copies/L) and Active Case Number for ALLA Sewershed',
#        col="Legend") + # legend title
#   scale_color_manual(labels = c("Active case number median", "RNA level"), 
#                      values = c("red", "blue")) + # edit legend labels 
#   geom_ribbon(aes(ymin=activecase_l, ymax=activecase_h), linetype=2, alpha=0.1)

# Try another plot
ggplot(alla_est, aes(x=day, y=activecase_m)) +
  geom_line() +
  labs(title= "Active Case Number for ALLA Sewershed") +
  geom_ribbon(aes(ymin=activecase_l, ymax=activecase_h), linetype=1, alpha=0.1, fill = "gray10", colour="gray")

predict(table[[1]]$coef1, data$gene)

# Use epyestim
library(EpiEstim)
res_parametric_si <- estimate_R(alla_est$activecase_m, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 2.6, 
                                  std_si = 1.5))
)

head(res_parametric_si$R)

# Plot results
plot(res_parametric_si, legend = FALSE)

# Line chart showing relationship between R estimate and RNA level
# ggplot(data, aes(x=mean_gene, y=meanr)) +
#   geom_line()

# Scatterplot
# ggplot(data, aes(x=mean_gene, y=meanr)) +
#   geom_point()


### Deconvolution of wastewater data by JS Huisman et al. ###
###### Deconvolve #####

deconv_result <- data.frame()
result <- data.frame()
for (incidence_var_i in c('weighted_avg_sars_cov2')){
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

# Plot number of cases
ggplot() +
  geom_point(data=deconv_result,
              aes(x=date, y=value, color='green')) +
  labs(title="Incidence Estimated by Deconvolution Method (SW)") +
  ylab("Incidence") +
  theme(legend.position = "none")

