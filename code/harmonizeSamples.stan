data {
  int<lower=0> N;
  vector[N] y;
  int<lower=0> Nbottles;  //N unique bottle samples
  int<lower=0> Ninstances; //N unique sampling instances in space and time
  int<lower=0> Nbetas; //N beta params to be estimated
  int<lower=0> bottle_idx[N]; //unique bottle samples
  int<lower=0> instance_idx[Nbottles]; //unique sampling instances
  vector[Ninstances] logVol; //log volume sampled for each unique instance
  matrix[Ninstances, Nbetas] X; //design matrix; sampling instances (rows) and parameters (cols)
}

parameters {
  vector[Nbottles] mu1; //bottle-level mean (= mean of tech reps within a bottle)
  // real<lower=0> sigma1; //tech-rep-level SD
  real<lower=0> sigma2; //bottle-level SD
  vector[Nbetas] B; //vector of intercepts and treatment effects
  real alpha; //scale param on relationship between bottle-level mean and SD
  real theta; //slope param on relationship between bottle-level mean and SD
}

transformed parameters{
  vector<lower=0>[Nbottles] sigma1; //tech-rep-level SD
  vector[Ninstances] mu2; // mean across bottles for the same sampling instance + same treatment

  
  sigma1 = inv_logit(alpha + theta*mu1); //assumes max of SD = 1
  // sigma1 = alpha * exp(-1*theta*mu1);
  

  mu2 = X*B + logVol; //main linear regression, in matrix format; sampling volume as offset

}

model {
  for (n in 1:N){
    y[n] ~ normal(mu1[bottle_idx[n]], sigma1[bottle_idx[n]]);
  }
  
  for (i in 1:Nbottles){
    mu1[i] ~ normal(mu2[instance_idx[i]], sigma2);
  }
  
  //priors
  mu1 ~ normal(0,10);
  // sigma1 ~ gamma(1,1);
  sigma2 ~ gamma(1,1);
  B ~ normal(0,5);
  alpha ~ normal(0,5);
  theta ~ normal(0,2);
}



