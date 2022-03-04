library(tidyverse)

# ---- Read data ---- 

E1_data <- read.csv('data/Exp1.csv')
E2_data <- read.csv('data/Exp2.csv')
E3_data <- read.csv('data/Exp3.csv')

# ---- Calculate distractor interference ----

E1_sdata <- group_by(E1_data, dist_cond = if_else(as_pos>0, if_else(as_updown==1, "Freq", "Rare"), "Absent"), 
         group, dtfreq, exptype, ratio, sub) %>% 
  mutate(distonborder = ifelse(dist_cond=="Absent", "No_dist", borderpos)) %>%
  filter(correct==1, rt > 0.15, rt < 3, tno > 120, distonborder!="Border") %>%
  summarize(rt=mean(rt, na.rm=TRUE)) %>% spread(dist_cond, rt) %>% 
  mutate(pc=Rare-Freq, npc=pc/Absent, mint=(Freq+Rare)/2-Absent, nmint=mint/Absent,
         fint=Freq-Absent, nfint=fint/Absent, rint=Rare-Absent, nrint=rint/Absent)

E1_sdata$exp <- "Exp. 1"
names(E1_sdata)[5] <- "participant"

E2_sdata <- group_by(E2_data, dist_cond, 
         group, dtfreq, exptype, ratio, participant) %>%
  mutate(distonborder = ifelse(dist_cond=="Absent", "No_dist", borderpos)) %>%
  filter(response_key.corr>0, response_key.rt > 0.15, response_key.rt < 3,
         trials.thisTrialN>119, distonborder!="Border") %>%
  summarize(rt=mean(response_key.rt, na.rm=TRUE)) %>% spread(dist_cond, rt) %>% 
  mutate(pc=Rare-Freq, npc=pc/Absent, mint=(Freq+Rare)/2-Absent, nmint=mint/Absent,
         fint=Freq-Absent, nfint=fint/Absent, rint=Rare-Absent, nrint=rint/Absent)

E2_sdata$exp <- "Exp. 2"

E3_sdata <- group_by(E3_data, dist_cond, 
                     group, dtfreq, exptype, ratio, participant) %>%
  mutate(distonborder = ifelse(dist_cond=="Absent", "No_dist", borderpos)) %>%
  filter(response_key.corr>0, response_key.rt > 0.15, response_key.rt < 3, 
         trials.thisTrialN>119, distonborder!="Border") %>%
  summarize(rt=mean(response_key.rt, na.rm=TRUE)) %>% spread(dist_cond, rt) %>% 
  mutate(pc=Rare-Freq, npc=pc/Absent, mint=(Freq+Rare)/2-Absent, nmint=mint/Absent,
         fint=Freq-Absent, nfint=fint/Absent, rint=Rare-Absent, nrint=rint/Absent)

E3_sdata$exp <- "Exp. 3"

comb_sdata <- rbind(E1_sdata, E2_sdata, E3_sdata)

comb_ssdata <- comb_sdata %>% summarize(mfint=mean(fint)*1000, 
                                        sfint=sd(fint/sqrt(n()-1))*1000,
                                        mrint=mean(rint)*1000, 
                                        srint=sd(rint/sqrt(n()-1))*1000)

# ---- Reorganize data ----
# Reorganize data into a single column with distractor inteference
# in both the frequent and the rare region

sr_data <- comb_ssdata  %>%
  pivot_longer(cols = c(mfint, mrint, sfint, srint), 
               names_to =  c(".value", "region"), names_pattern="(.)(.)")
names(sr_data)[6:7] <- c("mint", "seint")

# Calculate the local region distractor frequency
sr_data <- sr_data  %>% 
  mutate(drfreq=ifelse(region=="f", dtfreq*as.numeric(substr(ratio,1,2))/100,
                       dtfreq*as.numeric(substr(ratio,4,5))/100))
sr_data$drfreq <- round(sr_data$drfreq*200)/200 # Avoid rounding error problem
sr_data$dtfreq <- factor(sr_data$dtfreq)

# ---- Model fitting ----

# linear model 

fit_fun_linear <- function(x, C) {
  y=C[1]-C[2]*x
  return(y)
}

linear_fit <- function(C) {
  return(sum((sr_data$mint - fit_fun_linear(sr_data$drfreq,C))^2))
}

linear_fit_par <- optim(c(300, 300), linear_fit)


# Linear information model

fit_fun_log <- function(x, C) {
  y=C[1]-C[2]*log(x)
  return(y)
}

log_fit <- function(C) {
  return(sum((sr_data$mint - fit_fun_log(sr_data$drfreq,C))^2))
}

log_fit_par <- optim(c(100,100), log_fit)


# Exponential model

fit_fun_exp <- function(x, C) {
  y=C[1]*exp(-C[2]*x)
  return(y)
}

exp_fit <- function(C) {
  return(sum((sr_data$mint - fit_fun_exp(sr_data$drfreq,C))^2))
}

exp_fit_par <- optim(c(200,1), exp_fit)


# Normalized exponential model

fit_fun_exp2 <- function(x, C) {
  y=C[1]*(C[3]*exp(-C[4]*x))/(C[2]+C[3]*exp(-C[4]*x))
  return(y)
}

exp2_fit <- function(C) {
  return(sum((sr_data$mint - fit_fun_exp2(sr_data$drfreq,C))^2))
}

exp2_fit_par <- optim(c(300,1,1,1), exp2_fit)

fit_fun_pow <- function(x, C) {
  y=C[1]*(x)^(C[2])
  return(y)
}

pow_fit <- function(C) {
  return(sum((sr_data$mint - fit_fun_pow(sr_data$drfreq,C))^2))
}

pow_fit_par <- optim(c(100,-0.1), pow_fit)


# Normalized information model

fit_fun_log2 <- function(x, C) {
  y=C[1]*(-C[3]*log(x))/(C[2]-C[3]*log(x))
  return(y)
}

log2_fit <- function(C) {
  return(sum((sr_data$mint - fit_fun_log2(sr_data$drfreq,C))^2))
}

log2_fit_par <- optim(c(300,1,1), log2_fit)


N <- nrow(sr_data)

# Calculate AIC
# Exponential has 2 parameters, same as linear
exp_delta_AIC <- N*log(exp_fit_par$value) - N*log(linear_fit_par$value)

# Normalized exponential has 4 parameters, two more than linear
exp_norm_delta_AIC <- N*log(exp2_fit_par$value) - N*log(linear_fit_par$value) + 2*2

# Linear information has 2 parameters, same as linear
lin_inf_delta_AIC <- N*log(log_fit_par$value) - N*log(linear_fit_par$value) 

# Exponential information has 2 parameters, same as linear
pow_delta_AIC <- N*log(pow_fit_par$value) - N*log(linear_fit_par$value) 

# Normalized information has 3 parameters, 1 more than linear
inf_norm_delta_AIC <- N*log(log2_fit_par$value) - N*log(linear_fit_par$value) + 2


fit_fun <- function(x) {fit_fun_log(x,log_fit_par$par)}

mod_fig <- ggplot(sr_data, aes(x=drfreq, y=mint, color=dtfreq)) + geom_point() +
  geom_errorbar(aes(ymin=mint-seint, ymax=mint+seint)) + theme_bw() +
  labs(x="Local distractor frequency", y="Distractor interference (ms)", 
       color="Distractor prevalence") +  geom_function(fun=fit_fun, color="Black") + 
  theme(legend.position = "bottom") + scale_x_continuous(breaks=c(0.050, 0.100, 0.150, 0.200, 0.25, 0.300,
                                                                  0.350, 0.400, 0.450),
                                                         labels = c("5%", "10%", "15%", "20%", "25%", 
                                                                    "30%", "35%", "40%", "45%"))

ggsave("figures/mod.png", mod_fig, width=6, height=5)

# ---- Reorganize data - normalized  interference ----
# Reorganize data into a single column with distractor inteference
# in both the frequent and the rare region

comb_ssdata_norm <- comb_sdata %>% summarize(mfint=mean(nfint), 
                                        sfint=sd(nfint)/sqrt(n()-1),
                                        mrint=mean(nrint), 
                                        srint=sd(nrint)/sqrt(n()-1))
sr_data_norm <- comb_ssdata_norm %>%
  pivot_longer(cols = c(mfint, mrint, sfint, srint), 
               names_to =  c(".value", "region"), names_pattern="(.)(.)")
names(sr_data_norm)[6:7] <- c("mint", "seint")

# Calculate the local region distractor frequency
sr_data_norm <- sr_data_norm  %>% 
  mutate(drfreq=ifelse(region=="f", dtfreq*as.numeric(substr(ratio,1,2))/100,
                       dtfreq*as.numeric(substr(ratio,4,5))/100))
sr_data_norm$drfreq <- round(sr_data_norm$drfreq*200)/200 # Avoid rounding error problem
sr_data_norm$dtfreq <- factor(sr_data_norm$dtfreq)


# ---- Model fitting - normalized interference  ----


log_fit_norm <- function(C) {
  return(sum((sr_data_norm$mint - fit_fun_log(sr_data_norm$drfreq,C))^2))
}

log_fit_par_norm <- optim(c(1,1), log_fit_norm)


fit_fun_norm <- function(x) {fit_fun_log(x,log_fit_par_norm$par)}

mod_fig_norm  <- ggplot(sr_data_norm, aes(x=drfreq, y=mint, color=dtfreq)) + geom_point() +
  geom_errorbar(aes(ymin=mint-seint, ymax=mint+seint)) + theme_bw() +
  labs(x="Local distractor frequency", y="Normalized distractor interference", 
       color="Distractor prevalence") +  geom_function(fun=fit_fun_norm, color="Black") + 
  theme(legend.position = "bottom") + scale_x_continuous(breaks=c(0.050, 0.100, 0.150, 0.200, 0.25, 0.300,
                                                                  0.350, 0.400, 0.450),
                                                         labels = c("5%", "10%", "15%", "20%", "25%", 
                                                                    "30%", "35%", "40%", "45%"))

ggsave("figures/mod_norm.png", mod_fig_norm, width=6, height=5)

