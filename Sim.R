# ------------------------------------------
# title: Simulation of COVID-19 using SEIR model and relationship with
#         wastewater RNA level
# ------------------------------------------

# install.packages('writexl')
library(tidyverse)
library(data.table)
library(writexl)

population <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Population Over Time.csv")
population_estimate <- population %>% 
  group_by(Region) %>%
  summarise(max_pop = max(catchment_population_ons_mid_2019))

# Define empty list as output
table <- list()

# Set plot margins to avoid "figure margins too large" error
par(mar=c(4,4,4,2))

# Beta = 0.15 i.e. R0 = 1.5
table[[1]] <- MC.COVID19.wastewater(Sim=100,
                                    Tm=153, 
                                    beta.s= 0.15, 
                                    gamma.e = 0.2,
                                    gamma.i = 0.1,
                                    p = 0.0001, 
                                    N = 100000,
                                    mu.V.max = 7.6,
                                    sd.V.max = 0.8,
                                    mu.V.20 = 3.5,
                                    sd.V.20 = 0.35,
                                    T.V.max = 5,
                                    Ts = 1.1,
                                    Temp = 16,
                                    mu.tau0 = 130,
                                    sd.tau0 = 25,
                                    mu.Q = 2.5,
                                    sd.Q = 0.15,
                                    G.mean = 128,
                                    G.sd = 10
)

# Beta = 0.2 i.e. R0 = 2
table[[2]] <- MC.COVID19.wastewater(Sim=100, 
                                    Tm=153, 
                                    beta.s= 0.2, 
                                    gamma.e = 0.2,
                                    gamma.i = 0.1,
                                    p = 0.001, 
                                    N = 100000,
                                    mu.V.max = 7.6,
                                    sd.V.max = 0.8,
                                    mu.V.20 = 3.5,
                                    sd.V.20 = 0.35,
                                    T.V.max = 5,
                                    Ts = 1.1,
                                    Temp = 16,
                                    mu.tau0 = 130,
                                    sd.tau0 = 25,
                                    mu.Q = 2.5,
                                    sd.Q = 0.15,
                                    G.mean = 128,
                                    G.sd = 10
)

### Very long run time
# Beta = 0.3 i.e. R0 = 3
table[[3]] <- MC.COVID19.wastewater(Sim=100, 
                                    Tm=153, 
                                    beta.s= 0.3, 
                                    gamma.e = 0.2,
                                    gamma.i = 0.1,
                                    p = 0.0001,  
                                    N = 100000,
                                    mu.V.max = 7.6,
                                    sd.V.max = 0.8,
                                    mu.V.20 = 3.5,
                                    sd.V.20 = 0.35,
                                    T.V.max = 5,
                                    Ts = 1.1,
                                    Temp = 16,
                                    mu.tau0 = 130,
                                    sd.tau0 = 25,
                                    mu.Q = 2.5,
                                    sd.Q = 0.15,
                                    G.mean = 128,
                                    G.sd = 10
)

# Extract values for 95% CI, 75%CI and median of case numbers estimated from RNA level
c_compartment_r1 <- table[[1]]$model$C
e_compartment_r1 <- table[[1]]$model$E
i_compartment_r1 <- table[[1]]$model$I

c_compartment_r2 <- table[[2]]$model$C
e_compartment_r2 <- table[[2]]$model$E
i_compartment_r2 <- table[[2]]$model$I

c_compartment_r3 <- table[[3]]$model$C
e_compartment_r3 <- table[[3]]$model$E
i_compartment_r3 <- table[[3]]$model$I

est_values_r1 <- as.data.frame(cbind(table[[1]]$x, table[[1]]$output1, table[[1]]$output2, 
                                     table[[1]]$output3, table[[1]]$output4, table[[1]]$output5))
names(est_values_r1) <- c("x", "y1", "y2", "y3", "y4", "y5")

est_values_r2 <- as.data.frame(cbind(table[[2]]$x, table[[2]]$output1, table[[2]]$output2, 
                                     table[[2]]$output3, table[[2]]$output4, table[[2]]$output5))
names(est_values_r2) <- c("x", "y1", "y2", "y3", "y4", "y5")

est_values_r3 <- as.data.frame(cbind(table[[3]]$x, table[[3]]$output1, table[[3]]$output2, 
                                     table[[3]]$output3, table[[3]]$output4, table[[3]]$output5))
names(est_values_r3) <- c("x", "y1", "y2", "y3", "y4", "y5")

# rna times and C compartment
incidence_var_r1 <- as.data.frame(cbind(c_compartment_r1, e_compartment_r1, i_compartment_r1))
incidence_var_r2 <- as.data.frame(cbind(c_compartment_r2, e_compartment_r2, i_compartment_r2))
incidence_var_r3 <- as.data.frame(cbind(c_compartment_r3, e_compartment_r3, i_compartment_r3))

# Flatten list of lists (table[[1]] to table[[3]])
# list_sep <- unlist(table, recursive = FALSE)

# Write each table into a separate excel file to fix
# mylists <- list('r1sim'=table[[1]][1:6], 
#                 'r2sim'=table[[2]][1:6], 
#                 'r3sim'=table[[3]][1:6])
# purrr::walk(names(mylists),
#             function(x){
#               writexl::write_xlsx(mylists[[x]],
#                         path = paste0(x, ".xlsx"))
#               })

# Pivot into xy table for plot
# reshape2 is deprecated; use pivot_longer instead of melt
mdf1 <- est_values_r1 %>% pivot_longer(y1:y5) %>%
  distinct()
mdf2 <- est_values_r2 %>% pivot_longer(y1:y5) %>%
  distinct()
mdf3 <- est_values_r3 %>% pivot_longer(y1:y5) %>%
  distinct()

# Plot the smoothed curves from estimated cases, not exactly the spline function
ggplot(data=mdf1, aes(x=x, y=value, group=name, color=name)) +
  geom_line() + # this line function smooths out all the estimated cases 
  xlim(0, 1e+13) + 
  ylim(0, 2500) +
  labs(title= "Estimated active case number and RNA level (N=1,000,000, R0=1.5)",
       col="Percentile") + # legend title
  scale_color_manual(labels = c("2.5%", "12.5%", "50%", "87.5%", "97.5%"), 
                     values = c("yellow", "orange", "red", "orange", "yellow")) + # edit legend labels
  xlab("Gene copies per day") +
  ylab("Number of cases")

ggplot(data=mdf2, aes(x=x, y=value, group=name, color=name)) +
  geom_line() + # this line function smooths out all the estimated cases 
  xlim(0, 9e+13) + 
  ylim(0, 20000) +
  labs(title= "Estimated active case number and RNA level (N=100,000, R0=2)",
       col="Percentile") + # legend title
  scale_color_manual(labels = c("2.5%", "12.5%", "50%", "87.5%", "97.5%"), 
                     values = c("yellow", "orange", "red", "orange", "yellow")) + # edit legend labels
  xlab("Gene copies per day") +
  ylab("Number of cases")

ggplot(data=mdf3, aes(x=x, y=value, group=name, color=name)) +
  geom_line() + # this line function smooths out all the estimated cases 
  xlim(0, 1.5e14) + 
  ylim(0, 40000) +
  labs(title= "Estimated active case number and RNA level (N=100,000, R0=3)",
       col="Percentile") + # legend title
  scale_color_manual(labels = c("2.5%", "12.5%", "50%", "87.5%", "97.5%"), 
                     values = c("yellow", "orange", "red", "orange", "yellow")) + # edit legend labels
  xlab("Gene copies per day") +
  ylab("Number of cases")

# Export results
write.csv(est_values_r1_pmod,"model_estimates_r1_pmod.csv", row.names=FALSE)
write.csv(incidence_var_r1_pmod,"incidence_variables_r1_pmod.csv", row.names=FALSE)

write.csv(est_values_r2,"model_estimates_r2.csv", row.names=FALSE)
write.csv(incidence_var_r2,"incidence_variables_r2.csv", row.names=FALSE)

write.csv(est_values_r3,"model_estimates_r3.csv", row.names=FALSE)
write.csv(incidence_var_r3,"incidence_variables_r3.csv", row.names=FALSE)

table[[1]]
