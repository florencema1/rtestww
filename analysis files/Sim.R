# ------------------------------------------
# title: Simulation of COVID-19 using SEIR model and relationship with
#         wastewater RNA level
# ------------------------------------------

library(data.table)
table <- MC.COVID19.wastewater(Sim=30,
                      Tm=90,
                      beta.s= 1,
                      gamma.e = 0.2,
                      gamma.i = 0.1,
                      p = 0.005,
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
est_values1 <- table$y %*% table$coef1
est_values2 <- table$y %*% table$coef2
est_values3 <- table$y %*% table$coef3
est_values4 <- table$y %*% table$coef4
est_values5 <- table$y %*% table$coef5
est_values <- as.data.frame(cbind(table$x, est_values1, est_values2, est_values3, est_values4, est_values5))
names(est_values) <- c("x", "y1", "y2", "y3", "y4", "y5")

table$mod

# Export results
write.csv(table,"WasteWater_SEIR.csv", row.names=FALSE)
write.csv(est_values,"model_estimates.csv", row.names=FALSE)

# Pivot into xy table for plot
# reshape2 is deprecated; use pivot instead of melt
mdf <- est_values %>% pivot_longer(y1:y5) %>%
  distinct()

# Plot the smoothed curves 
ggplot(data=mdf, aes(x=x, y=value, group=name, color=name)) +
  geom_line()
