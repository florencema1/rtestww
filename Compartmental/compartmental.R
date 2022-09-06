# ------------------------------------------
# title: Simulation of COVID-19 using modified SEIR model
#        from COVID19WastewaterModel by McMahan et. al
# https://doi.org/10.1016/S2542-5196(21)00230-8
# repo: https://github.com/scwatson812/COVID19WastewaterModel
# ------------------------------------------

library(tidyverse)
library(data.table)
library(writexl)
library(patchwork)

# Define empty list as output
table <- list()

# Set plot margins to avoid "figure margins too large" error
par(mar=c(4,4,4,2))

# Set themes for plots
theme_set(theme_minimal())

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
# Extract estimated RNA levels
rna_r1_sim <- data.frame(table[[1]]$rna) # 547 rows
rna_r2_sim <- data.frame(table[[2]]$rna) # 307 rows
rna_r3_sim <- data.frame(table[[3]]$rna) # 307 rows

# Add time lag columns
rna_r1_lag <- sapply(1:ncol(rna_r1_sim), function(i) lag(rna_r1_sim[[i]], n=7)) %>%
  c()
rna_r2_lag <- sapply(1:ncol(rna_r2_sim), function(i) lag(rna_r2_sim[[i]], n=7)) %>%
  c()
rna_r3_lag <- sapply(1:ncol(rna_r3_sim), function(i) lag(rna_r3_sim[[i]], n=7)) %>%
  c()

# Pivot long
rna_r1_sim <- unlist(rna_r1_sim, use.names = FALSE)
rna_r2_sim <- unlist(rna_r2_sim, use.names = FALSE)
rna_r3_sim <- unlist(rna_r3_sim, use.names = FALSE)

# Bind colunms, remove NA and convert to numeric
rna_r1_sim <- data.frame(cbind(rna_r1_sim, rna_r1_lag, rep("r1", length(rna_r1_sim)))) %>%
  drop_na()
rna_r1_sim[, 1:2] <- lapply(rna_r1_sim[, 1:2], FUN = function(y){as.numeric(y)})
rna_r2_sim <- data.frame(cbind(rna_r2_sim, rna_r2_lag, rep("r2", length(rna_r2_sim)))) %>%
  drop_na()
rna_r2_sim[, 1:2] <- lapply(rna_r2_sim[, 1:2], FUN = function(y){as.numeric(y)})
rna_r3_sim <- data.frame(cbind(rna_r3_sim, rna_r3_lag, rep("r3", length(rna_r3_sim)))) %>%
  drop_na()
rna_r3_sim[, 1:2] <- lapply(rna_r3_sim[, 1:2], FUN = function(y){as.numeric(y)})

### Make a dataframe of values from real wastewater data
weighted <- read.csv("./data/all_var (before interpolation).csv")

# Normalise population to 100k to match with simulations
weighted$norm_gc <- weighted$weighted_avg_sars_cov2*400*100000

# Add column for t-7 days
rna_pre <- weighted %>%
  group_by(region) %>%
  mutate(lag_gc = lag(norm_gc, n = 7)) %>%
  select(norm_gc, lag_gc, region) 

# Change column names for joining
colnames(rna_r1_sim)=colnames(rna_r2_sim)=colnames(rna_r3_sim)=colnames(rna_pre) <- 
  c("rna_t", "rna_t_minus_7", "group")

# Plot
r1 <- ggplot(data=rna_r1_sim, aes(x=rna_t, y=rna_t_minus_7)) +
  geom_point(color='cadetblue2', size=1) +
  labs(title= "\u03b2 = 0.15") +
  xlab(NULL) +
  ylab("RNA Level on Day t-7 (gc/day)") +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(limits=c(0,2e14),
                     breaks=seq(0,2e14,1e14)) +
  scale_y_continuous(limits=c(0,2e14))

r2 <- ggplot(data=rna_r2_sim, aes(x=rna_t, y=rna_t_minus_7)) +
  geom_point(color='cadetblue3', size=1) +
  labs(title= "\u03b2 = 0.2") +
  xlab("RNA Level on Day t (gc/day)") +
  ylab(NULL) +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(limits=c(0,2e14),
                     breaks=seq(0,2e14,1e14)) +
  scale_y_continuous(limits=c(0,2e14))

r3 <- ggplot(data=rna_r3_sim, aes(x=rna_t, y=rna_t_minus_7)) +
  geom_point(color='cadetblue4', size=1) +
  labs(title= "\u03b2 = 0.3") +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(limits=c(0,2e14),
                     breaks=seq(0,2e14,1e14)) +
  scale_y_continuous(limits=c(0,2e14))

# Show three graphs horizontally
r1 + r2 + r3
ggsave("plots/Betas combined.png", width = 8, height = 3, units = "in")

# Plot real data
ggplot(data=rna_pre, aes(x=rna_t, y=rna_t_minus_7, color=group)) +
  geom_point() +
  labs(title= "SARS-CoV-2 RNA in Wastewater - Regional Data") +
  xlab("RNA Level on Day t (gc/day)") +
  ylab("RNA Level on Day t-7 (gc/day)") 

# Join all dataframes
sim_plot <- rbind(rna_r1_sim, rna_r2_sim, rna_r3_sim) %>%
  drop_na() 

g <- ggplot() +
  geom_point(data=sim_plot, aes(x=rna_t, y=rna_t_minus_7, alpha=0.2,
                                shape=group)) +
  geom_point(data = rna_pre, aes(x=rna_t, y=rna_t_minus_7, color=group, alpha=1)) +
  labs(title= "SARS-CoV-2 RNA in Wastewater - Simulations and Regional Data") +
  xlab("RNA Level on Day t (gc/day)") +
  ylab("RNA Level on Day t-7 (gc/day)") +
  labs(colour="Regions", shape="\u03b2 for Simulation") +
  scale_shape(labels = c("0.15","0.2","0.3")) +
  coord_cartesian(ylim = c(0, 2e13), xlim = c(0, 2e13)) 
g + guides(alpha='none') 

# Export plot and results
ggsave("plots/Compartmental combined.png", width = 8, height = 4.5, units = "in")

# Export simulation results
write_csv(table, "SEIR simulations.csv")
write_csv(sim_plot, "SEIR plot combined.csv")
