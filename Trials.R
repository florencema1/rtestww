# ----------------------------------------------------------------
# Candidate number: 211392
# ----------------------------------------------------------------
# title: SEIR model outputs
# ----------------------------------------------------------------
library(tidyverse)
library(linelist)
library(EpiEstim)

# Set ggplot theme
theme_set(theme_classic())


# Read in RNA level
data <- read.csv("/Linked Data.csv")
data$X <- NULL

real_case <- read.csv("/Users/FM/Public/HDS/Summer Project/local files/region_2022-07-01.csv")
real_r <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/r_val_cleaned.csv")
real_r$X <- NULL

# Read in estimated case numbers from COVID19WasteWater SEIR model
model_case_r1 <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/model_estimates_r1.csv")
model_factor_r1 <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/incidence_variables_r1.csv")

model_case_r2 <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/model_estimates_r2.csv")
model_factor_r2 <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/incidence_variables_r2.csv")

model_case_r3 <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/model_estimates_r3.csv")
model_factor_r3 <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/incidence_variables_r3.csv")

# Transform into date format
data$date <- lubridate::ymd(data$date)
real_case$date <- lubridate::ymd(real_case$date)
real_r$date <- lubridate::ymd(real_r$date)
max(data$date)
min(data$date) 

# Remove incidence number from non-study period
real_case <- with(real_case, real_case[(date >= "2020-09-11" & date <= "2021-02-10"), ])
real_r <- with(real_r, real_r[(date >= "2020-09-11" & date <= "2021-02-10"), ])

# Insert index number as 'day number' for joining tables later
# df <- tibble::rowid_to_column(data, "day")
# case <- tibble::rowid_to_column(case, "day")

data$active1 <- predict(loess(model_case_r1$y3~model_case_r1$x), data$gene_copies_per_day)
data$active2 <- predict(loess(model_case_r2$y3~model_case_r2$x), data$gene_copies_per_day)
data$active3 <- predict(loess(model_case_r3$y3~model_case_r3$x), data$gene_copies_per_day)

# Plot SW
data_sw <- data %>%
  filter(region=='south west') %>%
  select(date, active1, active2, active3)

plot_active <- data_sw %>% pivot_longer(active1:active3)

ggplot(plot_active, aes(x=date, y=value, group=name, color=name)) +
  geom_line() +
  labs(title= "South West Active Case Number (N = 100,000)",
       col="Percentile") + # legend title
  scale_color_manual(labels = c("beta 0.15", "beta 0.2", "beta 0.3"), 
                     values = c("red", "green", "blue")) + # edit legend labels

# Estimate case number from RNA level median
plot(table[[2]]$x, table[[2]]$coef3)


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


# Plot real incidence for comparison
op <- par(mfrow=c(3, 3))
by(real_case, real_case$areaName, function(b) plot(newCasesBySpecimenDate ~ date, b, 
                                                   ylim=c(0, max(real_case$newCasesBySpecimenDate)), 
                                                   main=unique(b$areaName)))
par(op)

op <- par(mfrow=c(3, 3))
by(real_r, real_r$variable, function(b) plot(value ~ date, b, 
                                             ylim=c(0.6, max(real_r$value)), 
                                             main=unique(b$variable)))
par(op)


# Incidence of SW only
real_case_sw <- real_case[real_case$areaName=='South West',]
ggplot() +
  geom_line(data=real_case_sw,
            aes(x=date, y=newCasesBySpecimenDate)) +
  labs(title="Incidence by Specimen Date (SW)") +
  ylab("Incidence") 

# Gene
# Plot real incidence for comparison
op <- par(mfrow=c(3, 3))
by(real_case, real_case$areaName, function(b) plot(newCasesBySpecimenDate ~ date, b, 
                                                   ylim=c(0, max(real_case$newCasesBySpecimenDate)), 
                                                   main=unique(b$areaName)))
par(op)

op <- par(mfrow=c(3, 3))
by(real_r, real_r$variable, function(b) plot(value ~ date, b, 
                                             ylim=c(0.6, max(real_r$value)), 
                                             main=unique(b$variable)))
par(op)


# Incidence of SW only
real_case_sw <- real_case[real_case$areaName=='South West',]
ggplot() +
  geom_line(data=real_case_sw,
            aes(x=date, y=newCasesBySpecimenDate)) +
  labs(title="Incidence by Specimen Date (SW)") +
  ylab("Incidence") 

# Gene
