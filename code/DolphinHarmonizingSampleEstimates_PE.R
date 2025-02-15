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
  filter(preservation != "later") %>% # remove samples preserved in later because protocol was not quite right 
  filter(vol_filtered == 3) %>% # for now while troubleshooting model, only include 3L filtered 
  filter(experiment == "p_e") %>% # for now while troubleshooting model, only include 3L filtered 
  unite(experiment, filter_size, preservation, extraction, col = "treatment", remove = FALSE) %>% 
  mutate(logDNA = log(Tt_total_conc),
         logVol = log(vol_filtered),
         pore_size_um = as.factor(filter_size)) %>% 
  mutate(preservation = case_when(preservation == "longmires" ~ "z_longmires",
                                  TRUE ~ preservation))

# create additional indices
#distinguish sampling instances (which might have many different treatments but share a common real concentration) vs. unique treatments
e$reality_idx <- match(e$experiment, unique(e$experiment))
e$bottle_idx <- match(e$sample, unique(e$sample))
e$treatment_idx <- match(e$treatment, unique(e$treatment))

#arrange unique treatments separate from replicated observations, for model matrix
data_small <- e %>% 
  select(reality_idx, treatment_idx, preservation, extraction) %>% 
  distinct() 


#create model matrix, w a different intercept for each shared reality/environment sampled
# modMat <- model.matrix(~ 0 + reality_idx + pore_size_um + preservation + extraction, data = data_small)
modMat <- model.matrix(~ preservation*extraction, data = data_small)


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
# plot(mod1, par = c("alpha", "theta"))
plot(mod1, par = c("sigma1", "sigma2")) #SD among tech replicates, SD among bottle replicates

summary(mod1, par = c("sigma1", "sigma2"))$summary

betasOut <- 
  data.frame(param = colnames(modMat),
             mean_est = summary(mod1, par = "B")$summary[,1],
             sd = summary(mod1, par = "B")$summary[,3],
             p.025 = summary(mod1, par = "B")$summary[,4],
             p.975 = summary(mod1, par = "B")$summary[,8],
             n_eff = summary(mod1, par = "B")$summary[,9],
             Rhat = summary(mod1, par = "B")$summary[,10])

## So here, I am comparing to 
# -80, PCI

#mean estimates for each instance*treatment
meanEst <- (modMat %*% betasOut$mean_est)
meanObs <- e %>% 
  group_by(treatment_idx) %>% 
  summarise(m = mean(logDNA)) %>% 
  pull(m)

#all good; model means estimate the observations well
plot(meanObs~meanEst)
abline(0,1)

e %>% 
  filter(reality_idx == 1) %>% 
  pull(logDNA)

#for any given treatment, easily get the mean estimate like this:
sum(modMat[which(data_small$treatment_idx == 2),] * betasOut$mean_est)

#full posterior of the betas
Bpost <- extract(mod1, par = "B")$B
dim(Bpost) #Niterations in rows, outcomes in columns

#sample full posteriors for a given instance
hist((Bpost %*% modMat[which(data_small$treatment_idx == 2),]))

#do this for all instances and store results
ModelEstimates <- as.data.frame(matrix(NA, nrow = length(unique(e$treatment_idx)), ncol = 7))
for (i in 1:nrow(ModelEstimates)){
  ModelEstimates[i,1] <- i
  
  tmp <- (Bpost %*% modMat[which(data_small$treatment_idx == i),]) 
  
  ModelEstimates[i,2:7] <- c(mean(tmp), quantile(tmp, probs = c(.025, .25, .5, .75, .975)))
}
colnames(ModelEstimates) <- c("treatment_idx", "mean_est", "p0.025","p0.25","p0.5","p0.75","p0.975")

modOut <- e %>% 
  left_join(data.frame(treatment_idx = 1:stan_data$Ntreatments,
                       treatment_mean = summary(mod1, par = "mu2")$summary[,1],
                       treatment_2.5 = summary(mod1, par = "mu2")$summary[,4],
                       treatment_97.5 = summary(mod1, par = "mu2")$summary[,8])) %>% 
  left_join(data.frame(bottle_idx = 1:stan_data$Nbottles,
                       bottle_mean = summary(mod1, par = "mu1")$summary[,1],
                       bottle_2.5 = summary(mod1, par = "mu1")$summary[,4],
                       bottle_97.5 = summary(mod1, par = "mu1")$summary[,8]))

mod_shape_color <- modOut %>% 
  select(treatment_idx, preservation,extraction) %>% 
  distinct()

# checking model performance 
p1 <- ModelEstimates %>% 
  left_join(mod_shape_color, by="treatment_idx") %>% 
  mutate(meanObs = e %>% group_by(treatment_idx) %>% summarise(m = mean(logDNA)) %>% pull(m)) %>% 
  mutate(sdObs = e %>% group_by(treatment_idx) %>% summarise(sd = sd(logDNA)) %>% pull(sd)) %>%
  mutate(seObs = sdObs / sqrt(9)) %>% 
  mutate(ciObs = 1.96*seObs) %>% 
  ggplot(aes(x = meanObs, y = mean_est,color=factor(preservation), shape=factor(extraction))) +
  geom_point(size=3) +
  geom_segment(aes(x = meanObs, xend = meanObs, y = p0.025, yend = p0.975)) +
  geom_segment(aes(x = meanObs-ciObs, xend = meanObs+ciObs, y = mean_est, yend = mean_est)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  theme_bw() +
  labs(x="Observation Mean (log copies/uL)", y="Modeled Mean (log copies/uL)",color="Preservation Method", shape="Extraction Method") +
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"), labels = c("-80", "Desiccation","Shield","Longmires")) + 
  scale_shape_manual(values=c(15,16,17),labels = c("PCI", "Qiagen","Zymo")) 

ggsave(here("figures","SUPFIG6_p_e_model_mean_obs.png"))

# Plot collapsing tech reps, bio reps, seeing difference between treatments 

treatment_idx_to_plot <- modOut %>% 
  select(treatment,treatment_idx, treatment_mean) %>% 
  distinct() %>% 
  arrange(desc(treatment_mean)) %>% 
  mutate(treatment_idx_plot = 1:10) %>% 
  select(treatment_idx, treatment_idx_plot)

p2 <- modOut %>% 
  left_join(treatment_idx_to_plot, by="treatment_idx") %>% 
  ggplot() + 
  # geom_point(aes(x=treatment_idx_plot, y=logDNA, color="Technical Replicate"), shape = 4, size=3, alpha=0.5) +
  # geom_point(aes(x=treatment_idx_plot, y=bottle_mean, color="Bottle Mean"), size=3, alpha=0.5) +
  # geom_segment(aes(x = treatment_idx_plot, xend = treatment_idx_plot, y = bottle_2.5, yend = bottle_97.5), color=wes_palette(n=4, name="Zissou1")[2]) +
  # geom_point(aes(x=treatment_idx_plot, y=treatment_mean, color="Treatment Mean"), size=3, alpha=0.5) +
  # geom_segment(aes(x = treatment_idx_plot, xend = treatment_idx_plot, y = treatment_2.5, yend = treatment_97.5), color=wes_palette(n=4, name="Zissou1")[4]) +
  geom_point(aes(x=treatment_idx_plot, y=logDNA, shape="Technical Replicate", color="Technical Replicate"), size=3, position=position_nudge(x = -0.2)) +
  geom_point(aes(x=treatment_idx_plot, y=bottle_mean, shape="Bottle Mean", color="Bottle Mean"), size=2, alpha=0.5) +
  geom_segment(aes(x = treatment_idx_plot, xend = treatment_idx_plot, y = bottle_2.5, yend = bottle_97.5), color=wes_palette(n=4, name="Zissou1")[2]) +
  geom_point(aes(x=treatment_idx_plot, y=treatment_mean, shape = "Treatment Mean", color="Treatment Mean"), size=4, alpha=0.5, position=position_nudge(x = 0.2)) +
  geom_segment(aes(x = treatment_idx_plot, xend = treatment_idx_plot, y = treatment_2.5, yend = treatment_97.5), color=wes_palette(n=4, name="Zissou1")[4], position=position_nudge(x = 0.2)) +
  labs(x="Treatment",y="Target DNA (log copies/uL)", color=" ") + 
  theme_bw() +
  scale_shape_manual(breaks= c("Technical Replicate", "Bottle Mean", "Treatment Mean"), values = c("Technical Replicate" = 4, "Bottle Mean" = 19, "Treatment Mean"=19), guide="none") +
  scale_color_manual(breaks= c("Technical Replicate", "Bottle Mean", "Treatment Mean"), values = c("Technical Replicate" = "darkgrey", "Bottle Mean" = wes_palette(n=4, name="Zissou1")[2], "Treatment Mean"=wes_palette(n=4, name="Zissou1")[4])) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("-80 / Q", "-80 / Z", "-80 / PCI","Shield / Z","Long / PCI", "Shield / Q","Shield / PCI","Desicc / Z","Desicc / Q","Desicc / PCI")) +
  theme(axis.text.x = element_text(angle = -90, vjust=0.5))
ggsave(here("figures","p_e_model_results.png"))

# plot_grid(p1, p2, labels = c('A', 'B'), nrow = 2, label_size = 12)
# ggsave(here("figures","p_e_model_results.png"), units="in", width=7, height=7)


###### CONVERTING
### 

