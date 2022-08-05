# ------------------------------------------
# title: Simulation of COVID-19 using SEIR model and relationship with
#         wastewater RNA level
# ------------------------------------------

# install.packages('writexl')
library(tidyverse)
library(data.table)
library(writexl)


# Define empty list as output
table <- list()

# Set plot margins to avoid "figure margins too large" error
par(mar=c(4,4,4,2))

# Beta = 0.15 i.e. R0 = 1.5
table[[1]] <- MC.COVID19.wastewater(Sim=100,
                                    Tm=547,
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
                                    Tm=307, 
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
                                    Tm=307, 
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

# Scatterplots of RNA over time for beta = 0.15/0.2/0.3
### Add wastewater RNA level on day (t-7)
# Convert RNA level to a single column first 
rna_r1_long <- unlist(data.frame(table[[1]]$rna))
rna_r2_long <- unlist(data.frame(table[[2]]$rna))
rna_r3_long <- unlist(data.frame(table[[3]]$rna))

# Add new columns
rna_r1_minus <- data.frame(cbind(rna_r1_long, lag(rep(rna_r1_long), n=7), 
                                 rep("r1", length(rna_r1_long)))) 
rna_r1_minus[, 1:2] <- lapply(rna_r1_minus[, 1:2], FUN = function(y){as.numeric(y)})

rna_r2_minus <- data.frame(cbind(rna_r2_long, lag(rep(rna_r2_long), n=7), 
                                 rep("r2", length(rna_r2_long)))) 
rna_r2_minus[, 1:2] <- lapply(rna_r2_minus[, 1:2], FUN = function(y){as.numeric(y)})

rna_r3_minus <- data.frame(cbind(rna_r3_long, lag(rep(rna_r3_long), n=7), 
                                 rep("r3", length(rna_r3_long)))) 
rna_r3_minus[, 1:2] <- lapply(rna_r3_minus[, 1:2], FUN = function(y){as.numeric(y)})

# Make a dataframe of values from real wastewater data
data <- read.csv("./Data/Linked Data.csv")
population <- read.csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Population Over Time.csv")
max_pop <- population %>%
  group_by(region) %>%
  summarise(max_pop = max(catchment_population_ons_mid_2019))

# Normalise to population 100k to match with simulations
rna_real_norm <- data %>%
  merge(max_pop, by='region')
rna_real_norm$gcpd_norm <- rna_real_norm$gene_copies_per_day/rna_real_norm$max_pop*100000

# Prepare data with RNA levels 7 days prior
rna_pre <- rna_real_norm %>%
  group_by(region) %>%
  mutate(lag_gc = lag(gcpd_norm, n = 7)) %>%
  select(gcpd_norm, lag_gc, region) %>%
  drop_na()

# Change column names for joining
colnames(rna_r1_minus)=colnames(rna_r2_minus)=colnames(rna_r3_minus)=colnames(rna_pre) <- 
  c("rna_t", "rna_t_minus_7", "group")



# Plot
ggplot(data=rna_r1_minus, aes(x=rna_t, y=rna_t_minus_7)) +
  geom_point(color='green') +
  labs(title= "Wastewater RNA Level of SARS-CoV-2 at Time(t) & Time(t-7), beta = 0.15") +
  xlab("RNA level at Time t") +
  ylab("RNA level at Time (t-7)")

ggplot(data=rna_r2_minus, aes(x=rna_t, y=rna_t_minus_7)) +
  geom_point(color='yellow') +
  labs(title= "Wastewater RNA Level of SARS-CoV-2 at Time(t) & Time(t-7), beta = 0.2") +
  xlab("RNA level at Time t") +
  ylab("RNA level at Time (t-7)")
  # geom_point(data = rna_pre, color = "red", shape=18) 

ggplot(data=rna_r3_minus, aes(x=rna_t, y=rna_t_minus_7)) +
  geom_point(color='blue') +
  labs(title= "Wastewater RNA Level of SARS-CoV-2 at Time(t) & Time(t-7), beta = 0.3") +
  xlab("RNA level at Time t") +
  ylab("RNA level at Time (t-7)")

ggplot(data=rna_pre, aes(x=rna_t, y=rna_t_minus_7, color=group)) +
  geom_point() +
  labs(title= "Wastewater RNA Level of SARS-CoV-2 at Time(t) & Time(t-7)") +
  xlab("RNA level at Time t") +
  ylab("RNA level at Time (t-7)")


# Join all dataframes
rna_minus <- rbind(rna_r1_minus, rna_r2_minus, rna_r3_minus, rna_pre) %>%
  drop_na()

ggplot(data=rna_minus, aes(x=rna_t, y=rna_t_minus_7, color=group)) +
  geom_point() +
  geom_point(data = rna_pre, aes(x=rna_t, y=rna_t_minus_7)) +
  labs(title= "Wastewater RNA Level of SARS-CoV-2 at Time(t) & Time(t-7)") +
  xlab("RNA level at Time t") +
  ylab("RNA level at Time (t-7)")
  # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
