# ------------------------------------------
# title: Simulation of COVID-19 using SEIR model and relationship with
#         wastewater RNA level
# ------------------------------------------
library(tidyverse)
library(data.table)

# Population has to be over a certain value?
pop <- c(50000,70000,80000)

# Define empty list as output
table <- list()

# for(j in 1:length(pop)) {
# table[[j]] <- MC.COVID19.wastewater(Sim=10,
#                       Tm=90,
#                       beta.s= 1,
#                       gamma.e = 0.2,
#                       gamma.i = 0.1,
#                       p = 0.05,
#                       N = pop[j],
#                       mu.V.max = 7.6,
#                       sd.V.max = 0.8,
#                       mu.V.20 = 3.5,
#                       sd.V.20 = 0.35,
#                       T.V.max = 5,
#                       Ts = 4,
#                       Temp = 18,
#                       mu.tau0 = 1,
#                       sd.tau0 = 0.1,
#                       mu.Q = 2.5,
#                       sd.Q = 0.15,
#                       G.mean = 128,
#                       G.sd = 10
# )
# }

# Single population for testing
table[[1]] <- MC.COVID19.wastewater(Sim=10,
                      Tm=90,
                      beta.s= 1,
                      gamma.e = 0.2,
                      gamma.i = 0.1,
                      p = 0.05,
                      N = 62500,
                      mu.V.max = 7.6,
                      sd.V.max = 0.8,
                      mu.V.20 = 3.5,
                      sd.V.20 = 0.35,
                      T.V.max = 5,
                      Ts = 4,
                      Temp = 18,
                      mu.tau0 = 1,
                      sd.tau0 = 0.1,
                      mu.Q = 2.5,
                      sd.Q = 0.15,
                      G.mean = 128,
                      G.sd = 10
)

# Extract values for 95% CI, 75%CI and median of case numbers estimated from RNA level
# est_values1 <- table[[1]]$y %*% table[[1]]$coef1
# est_values2 <- table[[1]]$y %*% table[[1]]$coef2
# est_values3 <- table[[1]]$y %*% table[[1]]$coef3
# est_values4 <- table[[1]]$y %*% table[[1]]$coef4
# est_values5 <- table[[1]]$y %*% table[[1]]$coef5
rna_time <- table[[1]]$rna
c_compartment <- table[[1]]$model$C
i_compartment <- table[[1]]$model$I
est_values <- as.data.frame(cbind(table[[1]]$x, est_values1, est_values2, est_values3, est_values4, est_values5))
names(est_values) <- c("x", "y1", "y2", "y3", "y4", "y5")

# rna times and C compartment
incidence_var <- as.data.frame(cbind(c_compartment, i_compartment))

# Flatten list of lists
list_sep <- flatten(table)

purrr::walk(names(table), 
            function(x){
              write.csv(table[[x]], 
                        path = paste0(x, ".xlsx"))
              })

# Pivot into xy table for plot
# reshape2 is deprecated; use pivot instead of melt
mdf <- est_values %>% pivot_longer(y1:y5) %>%
  distinct()

# Plot the smoothed curves *** from estimated cases not exactly the spline function
ggplot(data=mdf, aes(x=x, y=value, group=name, color=name)) +
  geom_line() + # this line function smooths out all the estimated cases 
  labs(title= "Estimated active case number and RNA level at N=60,000",
       col="Percentile") + # legend title
  scale_color_manual(labels = c("2.5%", "12.5%", "50%", "87.5%", "97.5%"), 
                     values = c("yellow", "orange", "red", "orange", "yellow")) + # edit legend labels
  xlab("Gene copies/L") +
  ylab("Number of cases")

# Test smoothed lines
plot <- data.frame(table[[1]]$x) %>%
  cbind(table[[1]]$coef1) %>%
  cbind(table[[1]]$coef2) %>%
  cbind(table[[1]]$coef3) %>%
  cbind(table[[1]]$coef4) %>%
  cbind(table[[1]]$coef5)
names(plot) <- c("x", "coef1", "coef2", "coef3", "coef4", "coef5")

plot_long <- plot%>%
  pivot_longer('coef1':'coef5')

# Plot
ggplot(plot_long, aes(x=x, y=value, group=name, color=name)) +
  geom_point() +
  scale_color_manual(labels = c("2.5%", "12.5%", "50%", "87.5%", "97.5%"), 
  values = c("red", "blue", "yellow", "green", "gray")) # edit legend labels 


# Export results
write.csv(est_values,"model_estimates.csv", row.names=FALSE)
write.csv(incidence_var,"incidence_variables.csv", row.names=FALSE)
write.csv(plot,"lines.csv", row.names=FALSE)

############################################################################
# Deconvolution

deconv_result <- data.frame()
result <- data.frame()
  new_deconv_result = deconvolveIncidence(clean_data, 
                                          incidence_var = gene,
                                          getCountParams('incubation'), 
                                          getCountParams('benefield'),
                                          smooth_param = TRUE, n_boot = 50) %>% #n_boot = 1000 in paper
    mutate(data_type = gene)
  
  new_result = getReBootstrap(new_deconv_result)
  new_result = new_result %>%
    mutate(data_type = gene)
  new_result['variable'] = gene
  
  deconv_result = bind_rows(deconv_result, new_deconv_result)
  result = bind_rows(result, new_result)
  
}
