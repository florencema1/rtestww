# ------------------------------------------
# title: Transforming RNA level into case number
# ------------------------------------------
library(tidyverse)

# Read in weekly RNA level for ALLA site
data <- read.csv("ALLA.csv")

# Read in R values of UK within study period
r_val <- read.csv("UK R.csv")

# Read in estimated case numbers from COVID19WasteWater SEIR model
case <- read.csv("WasteWater SEIR Est.csv")

# Insert index number i.e. week number to R values
r_val <- tibble::rowid_to_column(r_val, "week")

# Get the mean of R from upper and lower bounds
r_val$meanr <- (r_val$Upper.bound+r_val$Lower.bound)/2

# Combine RNA level and R estimate into one table
data <- merge(data, r_val, by = 'week')
data <- replace(data, is.na(data), 0)

# Estimate case number from RNA level of ALLA
bx <- bs(data$mean_gene, df=10)
est_med <- cbind(1,bx) %*% table$coef3 #coef 3 corresponds to median

# Create dataframe for plotting
alla_est <- as.data.frame(cbind(week = data$week, 
                              rna  = data$mean_gene)) %>%
          cbind(activecase=est_med) 
plotdf <- alla_est%>%
          pivot_longer('rna':'activecase')

write.csv(alla_est,"alla_est.csv", row.names=FALSE)

ggplot(plotdf, aes(x=week, y=value, group=name, color=name)) + 
  geom_line() +
  labs(title='RNA level (gene copies/L) and Active Case Number for ALLA Sewershed',
       col="Legend") + # legend title
  scale_color_manual(labels = c("Active case number", "RNA level"), 
                     values = c("#F8766D", "#619CFF")) # edit legend labels

# Use epyestim
library(EpiEstim)
res_parametric_si <- estimate_R(alla_est$activecase, 
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
