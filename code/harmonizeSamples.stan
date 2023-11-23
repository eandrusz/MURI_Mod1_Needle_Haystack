data {
  int<lower=0> N;
  vector[N] y;
  int<lower=0> Nbottles;  //N unique bottle samples
  int<lower=0> Ninstances; //N unique sampling instances in space and time
  int<lower=0> Nbetas; //N beta params to be estimated
  int<lower=0> bottle_idx[N];
  int<lower=0> instance_idx[Nbottles];
  vector[Ninstances] logVol; //log volume sampled for each unique instance
  matrix[Ninstances, Nbetas] X; //design matrix
}

parameters {
  vector[Nbottles] mu1;
  // vector[Ninstances] mu2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  vector[Nbetas] B;
}

transformed parameters{
  vector[Ninstances] mu2;

  mu2 = X*B + logVol;

}

model {
  for (n in 1:N){
    y[n] ~ normal(mu1[bottle_idx[n]], sigma1);
  }
  
  for (i in 1:Nbottles){
    mu1[i] ~ normal(mu2[instance_idx[i]], sigma2);
  }
  
  //priors
  
  mu1 ~ normal(0,10);
  // mu2 ~ normal(0,10);
  sigma1 ~ gamma(1,1);
  sigma2 ~ gamma(1,1);
  B[1] ~ normal(0, 10);
  B[2:Nbetas] ~ normal(0,2);
  
}



