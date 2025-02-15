#to harmonize sample estimates, given model of the effects of sampling volume, pore size, etc.
#RPK Nov 2023
#EAA May 2024 -- applying to volume/pore size and preservation/extraction for paper 
#EAA Jan 2025 -- editing, trying to remove volume and clean up 

library(here)
library(tidyverse)
select <- dplyr::select
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 4)
library(wesanderson)
library(cowplot)

a <- read_rds(here("data","experiments_merged.RDS"))

## Set up Stan model

e <- a %>% 
  #filter(vol_filtered == 3) %>% # for now while troubleshooting model, only include 3L filtered 
  filter(preservation == "longmires") %>% 
  filter(extraction == "PCI") %>% 
  # unite(vol_filtered, filter_size, preservation, extraction, col = "treatment", remove = FALSE) %>% 
  unite(experiment, vol_filtered, filter_size, col = "treatment", remove = FALSE) %>% 
  mutate(logDNA = log(Tt_total_conc),
         logVol = log(vol_filtered),
         pore_size_um = as.factor(filter_size)) 

# create additional indices
#distinguish sampling instances (which might have many different treatments but share a common real concentration) vs. unique treatments
e$reality_idx <- match(e$experiment, unique(e$experiment))
e$bottle_idx <- match(e$sample, unique(e$sample))
e$treatment_idx <- match(e$treatment, unique(e$treatment))

#arrange unique treatments separate from replicated observations, for model matrix
data_small <- e %>% 
  select(reality_idx, treatment_idx, pore_size_um, vol_filtered, logVol) %>% 
  distinct() %>% 
  mutate(reality_idx = as.factor(reality_idx)) #unique combinations of these variables


#track volume sampled
logVol <- data_small$logVol

#create model matrix, w a different intercept for each shared reality/environment sampled
# modMat <- model.matrix(~ 0 + reality_idx + pore_size_um + preservation + extraction, data = data_small)
modMat <- model.matrix(~ 0 + reality_idx + pore_size_um + logVol, data = data_small)


#make stan data object
stan_data <- list(
  N = nrow(e),
  y = e$logDNA,
  Nbottles = max(e$bottle_idx),
  Ntreatments = length(unique(e$treatment_idx)),
  # Nrealities = length(unique(e$reality_idx)),
  Nbetas = ncol(modMat),
  bottle_idx = e$bottle_idx,
  treatment_idx = e %>% 
    select(bottle_idx, treatment_idx) %>% 
    unique() %>% pull(treatment_idx),
  # reality_idx = e$reality_idx,
  #logVol = logVol,
  X = modMat #design matrix Ninstances x Nbetas
)

mod1 <- stan(file = here("code/harmonizeSamples_EA.stan"),
             data = stan_data,
             verbose = FALSE, chains = 3, thin = 1,
             warmup = 1000, iter = 2500,
             control = list(adapt_init_buffer = 175,
                            max_treedepth=12,
                            stepsize=0.01,
                            adapt_delta=0.7,
                            metric="diag_e"),
             # pars = stan_pars,
             refresh = 10,
             boost_lib = NULL 
) 

#shinystan::launch_shinystan(mod1)

plot(mod1, par = "B") #param fits
plot(mod1, par = "mu1") #bottle-level means
plot(mod1, par = "mu2") #treatment-level means
#plot(mod1, par = c("alpha", "theta"))
plot(mod1, par = c("sigma1", "sigma2")) #SD among tech replicates, SD among bottle replicates

summary(mod1, par = c("sigma1", "sigma2"))$summary

summary(mod1, par = "mu2")$summary[,1]

betasOut <- 
  data.frame(param = colnames(modMat),
             mean_est = summary(mod1, par = "B")$summary[,1],
             sd = summary(mod1, par = "B")$summary[,3],
             p.025 = summary(mod1, par = "B")$summary[,4],
             p.975 = summary(mod1, par = "B")$summary[,8],
             n_eff = summary(mod1, par = "B")$summary[,9],
             Rhat = summary(mod1, par = "B")$summary[,10])

realityEstimates <- betasOut[grep("reality", betasOut$param),]
treatmentEstimates <- betasOut[-grep("reality", betasOut$param),]

#mean estimates for each instance*treatment
meanEst <- (modMat %*% betasOut$mean_est) 
meanObs <- e %>% 
  #mutate(adjlogDNA = logDNA-logVol) %>% 
  group_by(treatment_idx) %>% 
  #summarise(m_adj = mean(adjlogDNA)) %>% 
  #pull(m_adj)
  summarise(m = mean(logDNA)) %>% 
  pull(m)

#all good; model means estimate the observations well
plot(meanObs~meanEst)
abline(0,1)

#full posterior of the betas
Bpost <- extract(mod1, par = "B")$B
dim(Bpost) #Niterations in rows, outcomes in columns

#sample full posteriors for a given instance
hist((Bpost %*% modMat[which(data_small$treatment_idx == 2),]))

modOut <- e %>% 
  left_join(data.frame(treatment_idx = 1:stan_data$Ntreatments,
                       treatment_mean = summary(mod1, par = "mu2")$summary[,1],
                       treatment_2.5 = summary(mod1, par = "mu2")$summary[,4],
                       treatment_97.5 = summary(mod1, par = "mu2")$summary[,8])) %>% 
  left_join(data.frame(bottle_idx = 1:stan_data$Nbottles,
                       bottle_mean = summary(mod1, par = "mu1")$summary[,1],
                       bottle_2.5 = summary(mod1, par = "mu1")$summary[,4],
                       bottle_97.5 = summary(mod1, par = "mu1")$summary[,8]))

#do this for all instances and store results
ModelEstimates <- as.data.frame(matrix(NA, nrow = length(unique(e$treatment_idx)), ncol = 7))
for (i in 1:nrow(ModelEstimates)){
  ModelEstimates[i,1] <- i
  
  tmp <- (Bpost %*% modMat[which(data_small$treatment_idx == i),]) 
  
  ModelEstimates[i,2:7] <- c(mean(tmp), quantile(tmp, probs = c(.025, .25, .5, .75, .975)))
}
colnames(ModelEstimates) <- c("treatment_idx", "mean_est", "p0.025","p0.25","p0.5","p0.75","p0.975")

mod_shape_color <- modOut %>% 
  select(treatment_idx, experiment, pore_size_um, vol_filtered) %>% 
  distinct()

#posterior predictive check, for the means
p1 <- ModelEstimates %>% 
  left_join(mod_shape_color, by="treatment_idx") %>% 
  mutate(meanObs = e %>% group_by(treatment_idx) %>% summarise(m = mean(logDNA)) %>% pull(m)) %>%
  mutate(sdObs = e %>% group_by(treatment_idx) %>% summarise(sd = sd(logDNA)) %>% pull(sd)) %>%
  mutate(seObs = sdObs / sqrt(9)) %>% 
  mutate(ciObs = 1.96*seObs) %>% 
  ggplot(aes(x = meanObs, y = mean_est, shape= as_factor(vol_filtered), color=as_factor(pore_size_um))) +
  geom_point(size=3) +
  geom_segment(aes(x = meanObs, xend = meanObs, y = p0.025, yend = p0.975)) +
  geom_segment(aes(x = meanObs-ciObs, xend = meanObs+ciObs, y = mean_est, yend = mean_est)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  annotate("text", label="Campaign 1",x=9.4,y=8.5, size=4, color="black") +
  annotate("text", label="Campaign 1",x=10.7,y=9.9, size=4, color="black") +
  annotate("text", label="Campaign 2",x=12.1,y=12.7, size=4, color="black") +
  theme_bw() +
  #scale_color_manual(breaks= c("1", "3"), values = c("1" = wes_palette(n=4, name="Chevalier1")[1], "3"=wes_palette(n=4, name="Chevalier1")[4])) +
  scale_color_manual(breaks= c("1", "5"), values = c("1" = wes_palette(n=4, name="Chevalier1")[1], "5"=wes_palette(n=4, name="Chevalier1")[4])) +
  labs(x="Observation Mean (log copies/uL)", y="Modeled Mean (log copies/uL)", color = "Filter Pore Size (um)",shape="Volume Filtered (L)")

ggsave(here("figures","SUPFIG5_v_ps_model_mean_obs.png"))

# Plot collapsing tech reps, bio reps, seeing difference between treatments 

treatment_idx_to_plot <- modOut %>% 
  select(treatment,treatment_idx, treatment_mean) %>% 
  distinct() %>% 
  arrange(desc(treatment_mean)) %>% 
  mutate(treatment_idx_plot = 1:5) %>% 
  select(treatment_idx, treatment_idx_plot)

p2 <- modOut %>% 
  left_join(treatment_idx_to_plot, by="treatment_idx") %>% 
  ggplot() + 
  geom_point(aes(x=treatment_idx_plot, y=logDNA, shape="Technical Replicate", color="Technical Replicate"), size=3, position=position_nudge(x = -0.2)) +
  geom_point(aes(x=treatment_idx_plot, y=bottle_mean, shape="Bottle Mean", color="Bottle Mean"), size=2, alpha=0.5) +
  geom_segment(aes(x = treatment_idx_plot, xend = treatment_idx_plot, y = bottle_2.5, yend = bottle_97.5), color=wes_palette(n=4, name="Zissou1")[2]) +
  geom_point(aes(x=treatment_idx_plot, y=treatment_mean, shape = "Treatment Mean", color="Treatment Mean"), size=4, alpha=0.5, position=position_nudge(x = 0.2)) +
  geom_segment(aes(x = treatment_idx_plot, xend = treatment_idx_plot, y = treatment_2.5, yend = treatment_97.5), color=wes_palette(n=4, name="Zissou1")[4], position=position_nudge(x = 0.2)) +
  labs(x="Treatment",y="Target DNA (log copies/uL)", color=" ") + 
  theme_bw() +
  annotate("text", label="Campaign 2",x=2,y=12.65, size=4, color="black") +
  geom_segment(aes(x=1.5,xend=5.5,y=10.8,yend=10.8)) +
  annotate("text", label="Campaign 1",x=3.5,y=11.1, size=4, color="black") +
  scale_shape_manual(breaks= c("Technical Replicate", "Bottle Mean", "Treatment Mean"), values = c("Technical Replicate" = 4, "Bottle Mean" = 19, "Treatment Mean"=19), guide="none") +
  scale_color_manual(breaks= c("Technical Replicate", "Bottle Mean", "Treatment Mean"), values = c("Technical Replicate" = "darkgrey", "Bottle Mean" = wes_palette(n=4, name="Zissou1")[2], "Treatment Mean"=wes_palette(n=4, name="Zissou1")[4])) +
  scale_x_continuous(breaks=c(1,2,3,4,5), labels=c("5 um, 3 L", "1 um, 3 L", "5 um, 3 L", "1 um, 1 L", "5 um, 1 L")) +
  theme(axis.text.x = element_text(angle = -90))

ggsave(here("figures","v_ps_model_results.png"))

# plot_grid(p1, p2, labels = c('A', 'B'), nrow = 2, label_size = 12)
# ggsave(here("figures","v_ps_model_results.png"), units="in", width=7, height=7)


###### CONVERTING
# so from mod Mat -- select B for either reality index 1 or 2
# baseline pore size is 1 um
# log vol has to be added... ? 

### FIRST, check that I know what I am doing -- recreate what we have 
c2_5um_3L = betasOut$mean_est[2] + betasOut$mean_est[3] + betasOut$mean_est[4]*1.1

c1_5um_3L = betasOut$mean_est[1] + betasOut$mean_est[3] + betasOut$mean_est[4]*1.1
c1_5um_1L = betasOut$mean_est[1] + betasOut$mean_est[3] + betasOut$mean_est[4]*0
c1_1um_3L = betasOut$mean_est[1] + betasOut$mean_est[3]*0 + betasOut$mean_est[4]*1.1
c1_1um_1L = betasOut$mean_est[1] + betasOut$mean_est[3]*0 + betasOut$mean_est[4]*0

##### estimates all combinations not done from campaign 2 
c2_5um_1L = betasOut$mean_est[2] + betasOut$mean_est[3] + betasOut$mean_est[4]*0
c2_1um_1L = betasOut$mean_est[2] + betasOut$mean_est[3]*0 + betasOut$mean_est[4]*0
c2_1um_3L = betasOut$mean_est[2] + betasOut$mean_est[3]*0 + betasOut$mean_est[4]*1.1

