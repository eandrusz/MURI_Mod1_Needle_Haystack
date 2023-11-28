#to harmonize sample estimates, given model of the effects of sampling volume, pore size, etc.

library(here)
library(tidyverse)
select <- dplyr::select
library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  rstan_options(threads_per_chain = 4)

a <- read.csv(here("data/All_Mod1_Molecular_Data - all_data.csv"))
a$dolphin_DNA_copy_uL[a$dolphin_DNA_copy_uL == 0]<- 0.1 #make all values nonzero, to handle logs. We can make this better later. 

# example: effect of extraction*preservation
#     b <- a %>% 
#       filter(task == "matrix") %>% 
#       select(dolphin_DNA_copy_uL, extraction, preservation, volume_L) %>% 
#       mutate(logDNA = log(dolphin_DNA_copy_uL)) %>% 
#       mutate(extraction = case_when(extraction == "zymo" ~ "_zymo",
#                                     extraction != "zymo" ~ extraction))
#     
#     b %>% 
#       ggplot(aes(x = extraction, y = logDNA)) +
#         geom_point() +
#         facet_grid(~preservation)
#     
#     lm1 <- lm(logDNA~preservation*extraction, data = b)
#     summary(lm1)
#     
#     b %>% 
#       mutate(pred = predict(lm1)) %>% 
#       ggplot(aes(x = extraction, y = logDNA)) +
#       geom_point() +
#       facet_grid(~preservation) +
#       geom_point(aes(x = extraction, y = pred), color = "red")
# 
# # example: filter volume + pore size
#     d <- a %>% 
#       filter(task == "filter_vol") %>% 
#       select(volume_L, pore_size_um, dolphin_DNA_copy_uL) %>% 
#       mutate(logDNA = log(dolphin_DNA_copy_uL)) %>% 
#       mutate(std_logDNA = (logDNA-mean(logDNA))/sd(logDNA))
#     
#     d %>% 
#       ggplot(aes(x = volume_L, y = logDNA, color = as.factor(pore_size_um))) +
#         geom_point()
#     
#     lm2 <- lm(logDNA~volume_L*pore_size_um, data = d)
#     summary(lm2)    
    

    

## Stan model
  

    
  e <- a %>% 
    # filter(extra <1) %>% 
    # filter(task != "ESP") %>% 
    # filter(task != "decay") %>%
    # filter(task == "matrix") %>% 
    # filter(task != "variability") %>%
    # filter(task != "transect") %>%
    # filter(task %in% c("variability", "transect")) %>% 
    filter(!str_detect(reference, "SR_ESP")) %>% 
    mutate(reference = str_replace_all(reference, "_[1-6]$", "")) %>% 
    unite(volume_L, pore_size_um, preservation, extraction, col = "treatment", remove = FALSE) %>% 
    select(lab_ID,
           reference,
          treatment,
          common_reality,
           bio_rep,
           tech_rep,
           volume_L,
           pore_size_um,
           preservation,
           extraction,
           dolphin_DNA_copy_uL) %>% 
    mutate(logDNA = log(dolphin_DNA_copy_uL),
           logVol = log(volume_L),
           pore_size_um = as.factor(pore_size_um)) %>% 
    select(-c(dolphin_DNA_copy_uL, volume_L))
    
    e$reality_idx <- match(e$common_reality, unique(e$common_reality))
    e$bottle_idx <- match(e$lab_ID, unique(e$lab_ID))
    e$treatment_idx <- match(e$treatment, unique(e$treatment))
    
    e <- e %>% 
      unite(reality_idx, treatment_idx, col = "instance", remove = FALSE) %>% 
      mutate(instance_idx = match(instance, unique(instance)))
    
    #NOTE create separate intercept for each unique site/date/time, because we expect these to have a common real concentration of eDNA
    #distinguish sampling instances (which might have many different treatments but share a common real concentration) vs. unique treatments
    
    #arrange unique treatments separate from replicated observations, for model matrix

    data_small <- e %>% 
      select(reality_idx, instance_idx, treatment_idx, pore_size_um, preservation, extraction, logVol) %>% 
      distinct() %>% 
      mutate(reality_idx = as.factor(reality_idx)) #unique combinations of these variables
      
    logVol <- data_small$logVol
    
    modMat <- model.matrix(~ 0 + reality_idx + pore_size_um + preservation*extraction, data = data_small)    # 
    
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
                 warmup = 500, iter = 1500,
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
    plot(mod1, par = c("sigma1", "sigma2")) #SD among tech replicates, SD among bottle replicates
    
    betasOut <- 
    data.frame(param = colnames(modMat),
    mean_est = summary(mod1, par = "B")$summary[,1])
    
    realityEstimates <- betasOut[grep("reality", betasOut$param),]
    treatmentEstimates <- betasOut[-grep("reality", betasOut$param),]
    
    #example: suppose we want to know what reality/environment 1 would look like if we had 5um filters, shield, zymo, and 3L sample
    #mean estimate is:
    realityEstimates[1,2] + #reality/environment intercept
      treatmentEstimates[1,2] + #5um filter
      treatmentEstimates[5,2] + #shield preservation
      treatmentEstimates[8,2] + #zymo extraction
      treatmentEstimates[20,2] + #shield/zymo interaction
      log(3) #3L sample
    # our observation of that reality was lower across the board, because we were using Longmire's and PCI
    e %>% 
      filter(reality_idx == 1)
    
    
    #to do this as a sample from the full posterior:
    sampleExample <- unlist(extract(mod1, par = "B[1]")) + #reality/environment intercept
      unlist(extract(mod1, par = "B[61]")) + #5um filter
               unlist(extract(mod1, par = "B[65]"))+ #shield preservation
                        unlist(extract(mod1, par = "B[68]"))+ #zymo extraction
                                 unlist(extract(mod1, par = "B[80]")) +
      log(3) #3L sample
    hist(sampleExample)
    
    #now do this for each reality/environment, and store summary:
    ProjectedEstimates <- as.data.frame(matrix(NA, nrow = length(unique(e$reality_idx)), ncol = 7))
      for (i in 1:nrow(ProjectedEstimates)){
        ProjectedEstimates[i,1] <- i
        
        tmp <- unlist(extract(mod1, par = paste0("B[", i,"]"))) + #reality/environment intercept
          unlist(extract(mod1, par = "B[61]")) + #5um filter
          unlist(extract(mod1, par = "B[65]"))+ #shield preservation
          unlist(extract(mod1, par = "B[68]"))+ #zymo extraction
          unlist(extract(mod1, par = "B[80]")) + #shield/zymo interaction
          log(3) #3L sample
        
        ProjectedEstimates[i,2:7] <- c(mean(tmp), quantile(tmp, probs = c(.025, .25, .5, .75, .975)))
      }
    colnames(ProjectedEstimates) <- c("reality_idx", "mean_est", "p0.025","p0.25","p0.5","p0.75","p0.975")
    
    ProjectedEstimates <- e %>% 
      select(reality_idx, reference) %>% 
      distinct() %>% 
      left_join(ProjectedEstimates) 
    
    #example of the ESP sample time series
    ProjectedEstimates %>% 
      filter(grepl("ESP_WCR", reference)) %>% 
      ggplot(aes(x = reality_idx, y = mean_est)) +
        geom_point() +
        geom_segment(aes(x = reality_idx, xend = reality_idx, y = p0.025, yend = p0.975))
    
    #note that here, I am giving the 95% CI for the MEAN estimate of what we would see for a given sample under given conditions. To add observation variability, we would sample from a normal distribution with SD = sqrt(sigma2 + sigma1) ...(bottle variability plus technical variability)
  
      
    
    
 
    
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
        ggplot(aes(x = bottle_idx, y = logDNA)) +
          geom_point() +
          geom_point(aes(x = bottle_idx, y = bottle_mean), color = "red") +
          geom_point(aes(x = bottle_idx, y = instance_mean), color = "purple")
        
          #facet_wrap(~instance)
    
    colnames(modMat)
    
    
    
    
