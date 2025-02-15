#to harmonize sample estimates, given model of the effects of sampling volume, pore size, etc.
#RPK Nov 2023
#EAA May 2024 -- applying to volume/pore size and preservation/extraction for paper 

library(here)
library(tidyverse)
select <- dplyr::select
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 4)
library(wesanderson)

a <- read_rds(here("data","experiments_merged.RDS"))

## Set up Stan model

e <- a %>% 
  filter(! is.na(Tt_total_conc)) %>% 
  filter(preservation != "later") %>% 
  unite(vol_filtered, filter_size, preservation, extraction, col = "treatment", remove = FALSE) %>% 
  mutate(logDNA = log(Tt_total_conc),
         logVol = log(vol_filtered),
         pore_size_um = as.factor(filter_size)) %>% 
  select(-c(Tt_total_conc, vol_filtered))

# create additional indices
#distinguish sampling instances (which might have many different treatments but share a common real concentration) vs. unique treatments
e$reality_idx <- match(e$experiment, unique(e$experiment))
e$bottle_idx <- match(e$sample, unique(e$sample))
e$treatment_idx <- match(e$treatment, unique(e$treatment))
e <- e %>% 
  unite(reality_idx, treatment_idx, col = "instance", remove = FALSE) %>% 
  mutate(instance_idx = match(instance, unique(instance)))

#arrange unique treatments separate from replicated observations, for model matrix
data_small <- e %>% 
  select(reality_idx, instance_idx, treatment_idx, pore_size_um, preservation, extraction, logVol) %>% 
  distinct() %>% 
  mutate(reality_idx = as.factor(reality_idx)) #unique combinations of these variables

#track volume sampled
logVol <- data_small$logVol

#create model matrix, w a different intercept for each shared reality/environment sampled
modMat <- model.matrix(~ 0 + reality_idx + pore_size_um + preservation*extraction, data = data_small)

#make stan data object
stan_data <- list(
  N = nrow(e),
  y = e$logDNA,
  Nbottles = max(e$bottle_idx),
  Ninstances = length(unique(e$instance)),
  Nbetas = ncol(modMat),
  bottle_idx = e$bottle_idx,
  instance_idx = e %>% 
    select(bottle_idx, instance_idx) %>% 
    unique() %>% pull(instance_idx),
  logVol = logVol,
  X = modMat #design matrix Ninstances x Nbetas
)

mod1 <- stan(file = here("code/harmonizeSamples.stan"),
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
plot(mod1, par = "mu2") #bottle-level means
plot(mod1, par = c("alpha", "theta"))
plot(mod1, par = c("sigma1", "sigma2")) #SD among tech replicates, SD among bottle replicates

summary(mod1, par = c("sigma1", "sigma2"))$summary

summary(mod1, par = "mu2")$summary[,1]
e %>% filter(instance_idx == 1)

betasOut <- 
  data.frame(param = colnames(modMat),
             mean_est = summary(mod1, par = "B")$summary[,1])

realityEstimates <- betasOut[grep("reality", betasOut$param),]
treatmentEstimates <- betasOut[-grep("reality", betasOut$param),]

#write_csv(treatmentEstimates, file = here("model_output/DNA_treatment_effects.csv"))

## So here, I am comparing to 
# 1 um, RNALater, PCI


#mean estimates for each instance*treatment
meanEst <- (modMat %*% betasOut$mean_est) + data_small$logVol
meanObs <- e %>% 
  #mutate(adjlogDNA = logDNA-logVol) %>% 
  group_by(instance_idx) %>% 
  #summarise(m_adj = mean(adjlogDNA)) %>% 
  #pull(m_adj)
  summarise(m = mean(logDNA)) %>% 
  pull(m)

#all good; model means estimate the observations well
plot(meanObs~meanEst)
abline(0,1)

# #example: suppose we want to know what reality/environment 1 would look like if we had 5um filters, shield, zymo, and 3L sample
# #mean estimate is:
# realityEstimates[1,2] + #reality/environment intercept
#   treatmentEstimates[1,2] + #5um filter
#   treatmentEstimates[5,2] + #shield preservation
#   treatmentEstimates[8,2] + #zymo extraction
#   treatmentEstimates[20,2] + #shield/zymo interaction
#   log(3) #3L sample
# # our observation of that reality was lower across the board, because we were using Longmire's and PCI

e %>% 
  filter(reality_idx == 1) %>% 
  pull(logDNA)

#for any given instance, easily get the mean estimate like this:
sum(modMat[which(data_small$instance_idx == 2),] * betasOut$mean_est) + data_small$logVol[2]



# #to do this as a sample from the full posterior:
# sampleExample <- unlist(extract(mod1, par = "B[1]")) + #reality/environment intercept
#   unlist(extract(mod1, par = "B[61]")) + #5um filter
#   unlist(extract(mod1, par = "B[65]"))+ #shield preservation
#   unlist(extract(mod1, par = "B[68]"))+ #zymo extraction
#   unlist(extract(mod1, par = "B[80]")) +
#   log(3) #3L sample
# hist(sampleExample)

#full posterior of the betas
Bpost <- extract(mod1, par = "B")$B
dim(Bpost) #Niterations in rows, outcomes in columns

#sample full posteriors for a given instance
hist((Bpost %*% modMat[which(data_small$instance_idx == 2),]) + data_small$logVol[2])

#do this for all instances and store results
ModelEstimates <- as.data.frame(matrix(NA, nrow = length(unique(e$instance_idx)), ncol = 7))
for (i in 1:nrow(ModelEstimates)){
  ModelEstimates[i,1] <- i
  
  tmp <- (Bpost %*% modMat[which(data_small$instance_idx == i),]) + data_small$logVol[i]
  
  ModelEstimates[i,2:7] <- c(mean(tmp), quantile(tmp, probs = c(.025, .25, .5, .75, .975)))
}
colnames(ModelEstimates) <- c("instance_idx", "mean_est", "p0.025","p0.25","p0.5","p0.75","p0.975")

#posterior predictive check, for the means
ModelEstimates %>% 
  mutate(meanObs = e %>% group_by(instance_idx) %>% summarise(m = mean(logDNA)) %>% pull(m)) %>% 
  ggplot(aes(x = meanObs, y = mean_est)) +
  geom_point() +
  geom_segment(aes(x = meanObs, xend = meanObs, y = p0.025, yend = p0.975)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  theme_bw() +
  labs(x="Observation Mean (log copies/uL)", y="Modeled Mean (log copies/uL)")

ggsave(here("figures","model_mean_obs.png"))

#now do this for each instance and project onto desired treatment, and store summary:
ProjectedEstimates <- as.data.frame(matrix(NA, nrow = length(unique(e$instance_idx)), ncol = 7))
for (i in 1:nrow(ProjectedEstimates)){
  ProjectedEstimates[i,1] <- i
  
  #get reality_idx for a given sampling index
  r <- as.integer(data_small$reality_idx[which(data_small$instance_idx == i)])
  
  #set desired treatment_idx to project estimates onto
  tr <- 6  #here, zymo and shield and 5um and 3 L
  
  #create vector of effect indicators
  effectVec <- rep(0, ncol(modMat))
  effectVec[c(which(colnames(modMat) == paste0("reality_idx",r)), #reality idx
              which(modMat[which(data_small$treatment_idx == 6),] == 1)[-1])] <- 1 #treatment indices
  
  tmp <- (Bpost %*% effectVec) + log(3) #project onto 3L samples
  
  ProjectedEstimates[i,2:7] <- c(mean(tmp), quantile(tmp, probs = c(.025, .25, .5, .75, .975)))
}
colnames(ProjectedEstimates) <- c("instance_idx", "mean_est", "p0.025","p0.25","p0.5","p0.75","p0.975")

ProjectedEstimates <- e %>% 
  # select(instance_idx, reality_idx, reference, logDNA) %>% 
  distinct() %>% 
  left_join(ProjectedEstimates) 
#write_csv(ProjectedEstimates, file = here("model_output/projected_logDNA_concentrations.csv"))


# #example of the ESP sample time series
# ProjectedEstimates %>% 
#   filter(grepl("ESP_WCR", reference)) %>% 
#   # arrange(reference) %>% 
#   ggplot(aes(x = reality_idx, y = mean_est)) +
#   geom_point() +
#   geom_segment(aes(x = reality_idx, xend = reality_idx, y = p0.025, yend = p0.975))
# 
# #note that here, I am giving the 95% CI for the MEAN estimate of what we would see for a given sample under given conditions. To add observation variability, we would sample from a normal distribution with SD = sqrt(sigma2 + sigma1) ...(bottle variability plus technical variability)






#in the model matrix, our reference condition for the intercepts is pore size 1um, Longmire's, and PCI extraction at 1L
colnames(modMat)
unique(data_small$extraction)
#we could change the reference condition to be our preferred combination: 3um, shield, zymo, and 3L


modOut <- e %>% 
  left_join(data.frame(instance_idx = 1:stan_data$Ninstances,
                       instance_mean = summary(mod1, par = "mu2")$summary[,1])) %>% 
  left_join(data.frame(bottle_idx = 1:stan_data$Nbottles,
                       bottle_mean = summary(mod1, par = "mu1")$summary[,1]))




modOut %>% 
  filter(reality_idx==2) %>% 
  select(bottle_idx, logDNA, bottle_mean, instance_mean) %>% 
  rename(observed = logDNA) %>% 
  pivot_longer(-bottle_idx, names_to = "Source") %>% 
  ggplot(aes(x = bottle_idx, y = value, color = Source)) +
  geom_point() +
  ylab("Log copies/L (Dolphin eDNA)") +
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))


modOut %>% 
  filter(reality_idx==2) %>% 
  #select(bottle_idx, logDNA, bottle_mean, instance_mean) %>% 
  #rename(observed = logDNA) %>% 
  #pivot_longer(-bottle_idx, names_to = "Source") %>% 
  ggplot(aes(x = preservation, y = logDNA)) +
  geom_boxplot() +
  geom_hline(yintercept=12.67542) +
  #geom_point(data=ModelEstimates, aes(x=preservation))
  facet_grid(~extraction, scales="free") + 
  ylab("Log copies/L (Dolphin eDNA)") +
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))

modOut %>% 
  #filter(reality_idx==1) %>% 
  filter(reality_idx==2) %>% 
  select(bottle_idx, logDNA, bottle_mean, instance_mean) %>% 
  rename(observed = logDNA) %>% 
  pivot_longer(-bottle_idx, names_to = "Source") %>% 
  ggplot(aes(x = Source, y = value)) +
  geom_point() +
  geom_hline(yintercept=12.67542) +
  #geom_hline(yintercept=10.08660) +
  #geom_point(data=ModelEstimates, aes(x=preservation)) +
  #facet_grid(~extraction, scales="free") + 
  ylab("Log copies/L (Dolphin eDNA)") +
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))



ProjectedEstimates %>% 
  filter(reality_idx == 2) %>% 
  ggplot(aes(x = tech_rep, y = logDNA)) +
    geom_point() +
    geom_point(aes(x = tech_rep, y = mean_est), color = "red") +
    geom_segment(aes(x = tech_rep, xend = tech_rep, y = p0.025, yend = p0.975), color = "red") +
    facet_grid(preservation~extraction) +
    scale_x_continuous(breaks=c(1,2,3)) +
    theme_bw() +
    labs(x="Technical Replicate",y="Modeled Mean (log copies/uL)")
#ggsave(here("figures","modelbytreatment_v_ps.png"))

ProjectedEstimates %>% 
  #filter(reality_idx == 2) %>% 
  ggplot(aes(x = reality_idx, y = logDNA, color=extraction)) +
  geom_point(alpha=0.3, size=3) +
  geom_point(aes(x = reality_idx, y = mean_est), color = "red", size=3) +
  geom_segment(aes(x = reality_idx, xend = reality_idx, y = p0.025, yend = p0.975), color = "red") +
  #facet_grid(preservation~extraction) +
  scale_x_continuous(breaks=c(1,2)) +
  theme_bw() +
  labs(x="Common Reality",y="Modeled Mean (log copies/uL)")
#ggsave(here("figures","model_byreality.png"), width=6, height=4, units="in")

ProjectedEstimates %>% 
  filter(reality_idx == 2) %>% 
  ggplot(aes(x = extraction, y = logDNA, color=preservation)) +
  geom_point(alpha=0.3, size=3) +
  geom_point(aes(x = extraction, y = mean_est), color = "black", size=3) +
  geom_segment(aes(x = extraction, xend = extraction, y = p0.025, yend = p0.975), color = "black") +
  #facet_grid(preservation~extraction) +
  #scale_x_continuous(breaks=c(1,2)) +
  theme_bw() +
  labs(x="Extraction Method",y="Modeled Mean (log copies/uL)", color="Preservation Method")
ggsave(here("figures","p_e_model_color.png"), width=6, height=4, units="in")

ProjectedEstimates %>% 
  filter(reality_idx == 1) %>% 
  ggplot(aes(x = logVol, y = logDNA, color=filter_size)) +
  geom_point(alpha=0.3, size=3) +
  geom_point(aes(x = logVol, y = mean_est), color = "black", size=3) +
  geom_segment(aes(x = logVol, xend = logVol, y = p0.025, yend = p0.975), color = "black") +
  #facet_grid(preservation~extraction) +
  #scale_x_continuous(breaks=c(1,2)) +
  theme_bw() +
  labs(x="Log Volume Sampled",y="Modeled Mean (log copies/uL)", color="Filter Pore Size (um)")
ggsave(here("figures","v_ps_model_color.png"), width=6, height=4, units="in")


tt<- modOut %>% 
  #filter(reality_idx==1) %>% 
  filter(reality_idx==2) %>% 
  select(bottle_idx, logDNA, bottle_mean, instance_mean) %>% 
  rename(observed = logDNA) %>% 
  pivot_longer(-bottle_idx, names_to = "Source") 

xx <- unique(modOut %>% select(bottle_idx, preservation, extraction))

yy <- tt %>% 
  left_join(xx) 

yy$Source <- factor(yy$Source, levels = c("observed", "bottle_mean", "instance_mean"))

ggplot(yy, aes(x = Source, y = value, color=factor(preservation), shape=factor(extraction))) +
  geom_point(size=3) +
  theme_bw() + 
  #geom_hline(yintercept= realityEstimates[2,2] + 3) +
  #geom_hline(yintercept=10.08660) +
  #geom_point(data=ModelEstimates, aes(x=preservation)) +
  #facet_grid(~extraction, scales="free") + 
  labs(x="", y="Log copies/L (Dolphin eDNA)", shape ="Extraction Method", color="Preservation Method") +
  scale_x_discrete(labels=c("All Replicates", "Technical Replicates Averaged", "Biological Replicates Averaged"),
                   guide = guide_axis(angle = -45)) + 
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))
ggsave(here("figures","p_e_tech_bio_instance.png"), width=8, height=6, units="in")


### for reality 1 

TT<- modOut %>% 
  filter(reality_idx==1) %>% 
  #filter(reality_idx==2) %>% 
  select(bottle_idx, logDNA, bottle_mean, instance_mean) %>% 
  rename(observed = logDNA) %>% 
  pivot_longer(-bottle_idx, names_to = "Source") 

XX <- unique(modOut %>% select(bottle_idx, pore_size_um, logVol))

YY <- TT %>% 
  left_join(XX) %>% 
  mutate(vol_filtered=exp(logVol)) %>% 
  mutate(adj_value = value-logVol) 

YY$Source <- factor(YY$Source, levels = c("observed", "bottle_mean", "instance_mean"))

ggplot(YY, aes(x = Source, y = adj_value, color=factor(vol_filtered), shape=factor(pore_size_um))) +
  geom_point(size=3) +
  theme_bw() + 
  #geom_hline(yintercept= realityEstimates[2,2] + 3) +
  #geom_hline(yintercept=10.08660) +
  #geom_point(data=ModelEstimates, aes(x=preservation)) +
  #facet_grid(~extraction, scales="free") + 
  labs(x="", y="Log copies/L (Dolphin eDNA)", shape ="Filter Pore Size (um)", color="Volume Filtered (L)") +
  scale_x_discrete(labels=c("All Replicates", "Technical Replicates Averaged", "Biological Replicates Averaged"),
                   guide = guide_axis(angle = -45)) + 
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1")) 
ggsave(here("figures","v_ps_tech_bio_instance.png"), width=8, height=6, units="in")


# modOut %>%
#   ggplot(aes(x = reality_idx, y = logDNA, color = as.factor(treatment_idx))) +
#   geom_point() +
#   geom_point(aes(x = reality_idx, y = bottle_mean), color = "black") +
#   geom_point(aes(x = reality_idx, y = instance_mean), color = "purple")
# # facet_wrap(~reality_idx)


## What we could do now is use ALL the data corrected for the different methods and see if we can say something about why they are still different concentations
##### number of dolphins, water temperature, etc?


# a %>% 
#   filter(extraction=="PCI", preservation=="longmires", vol_filtered==3, filter_size==5) %>% 
#   ggplot(aes(x=sample, y=Tt_total_conc)) + 
#   geom_point() +
#   facet_wrap(~experiment, scales="free_x")


