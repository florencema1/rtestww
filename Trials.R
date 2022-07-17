# ------------------------------------------
# title: Transforming RNA level into case number
# ------------------------------------------
library(tidyverse)
library(linelist)

# Set ggplot theme
theme_set(theme_bw())

# Read in weekly RNA level for ALLA site
data <- read.csv("/Users/FM/Public/HDS/Summer Project/rtestww/ALLA_day.csv")

# Read in R values of UK within study period
r_val <- read.csv("/Users/FM/Public/HDS/Summer Project/rtestww/UK R.csv")

# Get the mean of R from upper and lower bounds
r_val$meanr <- (r_val$Upper.bound+r_val$Lower.bound)/2

# Transform into date format
data$date <- lubridate::ymd(data$date)
r_val$Date <- lubridate::dmy(r_val$Date)

# Subset R values of date range of the RNA level file
subs_r <- r_val[r_val$Date >= min(data$date) & r_val$Date <= (max(data$date)+2),]

# Interpolate into daily values
data_int <- approx(x=data$date,y=data$gene, xout=seq(min(data$date), max(data$date), "days"))
list <- approx(x=subs_r$Date,y=subs_r$meanr, xout=seq(min(subs_r$Date), max(subs_r$Date), "days"))

# Trim r_val interpolated list and convert to dataframe
df <- lapply(list, function(x) {x <- x[-c(91:92)]})
df1 <- data.frame(df)
names(df1) <- c("date", "mean_r")

# Convert data interpolated into dataframe
df2 <- data.frame(data_int)
names(df2) <- c("date", "rna")

# Read in estimated case numbers from COVID19WasteWater SEIR model
case <- read.csv("/Users/FM/Public/HDS/Summer Project/rtestww/model_estimates.csv")
model_factor <- read.csv("/Users/FM/Public/HDS/Summer Project/incidence_variables.csv")

# Insert index number as 'day number' for joining tables later
df1 <- tibble::rowid_to_column(df1, "day")
df2 <- tibble::rowid_to_column(df2, "day")
case <- tibble::rowid_to_column(case, "day")

# Combine RNA level and R estimate into one table
combined <- merge(df1, df2, by = 'day')

# Replace all missing values with 0
# combined <- replace(combined, is.na(combined), 0)

# Estimate case number from RNA level of ALLA
bx <- bs(combined$rna, df=10)
est_med_med <- cbind(1,bx) %*% table[[1]]$coef3 #coef 3 corresponds to median
est_med_low <- cbind(1,bx) %*% table[[1]]$coef1 #coef 1, 5 refer to 95% CI; coef 2, 4 refer to 75% CI
est_med_high <- cbind(1,bx) %*% table[[1]]$coef5 

# Create dataframe for plotting
alla_est <- as.data.frame(cbind(day = combined$day, 
                              rna  = combined$rna)) %>%
  cbind(activecase_m=est_med_med) %>%
  cbind(activecase_l=est_med_low) %>%
  cbind(activecase_h=est_med_high) 

# Export est. values for ALLA 
write.csv(alla_est,"alla_est.csv", row.names=FALSE)

# Pivot dataframe for plotting multiple lines on the same graph
plotdf <- alla_est%>%
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
