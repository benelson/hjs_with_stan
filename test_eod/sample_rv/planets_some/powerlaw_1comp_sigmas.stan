data {
    int<lower=0> N; // number of planets w/o incl measurements
    int<lower=0> Ni; // number of planets w/ incl measurements

    real<lower=0.> per[N];
    real<lower=0.> rad[N];
    real<lower=0.> Mpl[N];
    
    real<lower=0.> per_sigma[N];
    real<lower=0.> rad_sigma[N];
    real<lower=0.> Mpl_sigma[N];

    real<lower=0.> peri[Ni];
    real<lower=0.> radi[Ni];
    real<lower=0.> Mpli[Ni];

    real<lower=0.> per_sigmai[Ni];
    real<lower=0.> rad_sigmai[Ni];
    real<lower=0.> Mpl_sigmai[Ni];

    real<lower=0.> ind[N];
    real<lower=0.> indi[Ni];
}
parameters {
    real gamma; // power law index
    real<lower=0.> xl; // lower truncation value of x
    real<upper=30.> xu; // upper truncation value of x
    real<lower=-1., upper=1.> cosi[N]; // cos inclinations of all planets

    real<lower=0.> per_true[N];
    real<lower=-3.,upper=1.> log_rad_true[N];
    real<lower=-3.,upper=2.> log_Mpl_true[N];

    real<lower=0.> peri_true[Ni];
    real<lower=-3., upper=1.> log_radi_true[Ni];
    real<lower=-3., upper=2.> log_Mpli_true[Ni];
}
transformed parameters{
    real rad_true[N];
    real Mpl_true[N];
    
    real radi_true[Ni];
    real Mpli_true[Ni];

    for (n in 1:N)
    {
	rad_true[n] <- 10. ^ log_rad_true[n];
      	Mpl_true[n] <- 10. ^ log_Mpl_true[n];
    }

    for (i in 1:Ni)
    {      
      	radi_true[i] <- 10. ^ log_radi_true[i];
    	Mpli_true[i] <- 10. ^ log_Mpli_true[i];
    }
}
model {
    real sum_log_x;
    real sini_onethird;

    real RJtoAU;
    real MJtoMsun;

    real x[N];
    real xi[Ni];

    RJtoAU <- 0.000477894503;
    MJtoMsun <- 0.000954265748;
    sum_log_x <- 0.0;
    
//    per_true ~ uniform(0.1,100.);
//    peri_true ~ uniform(0.1,100.);

    per ~ normal(per_true, per_sigma);
    rad ~ normal(rad_true, rad_sigma);
    Mpl ~ normal(Mpl_true, Mpl_sigma);

    peri ~ normal(peri_true, per_sigmai);
    radi ~ normal(radi_true, rad_sigmai);
    Mpli ~ normal(Mpli_true, Mpl_sigmai);
  
//    cosi ~ uniform(-1.,1.); // prior in cosi

    for (n in 1:N)
    {
        x[n] <- 0.462 * (per_true[n]/365.242)^(0.66667) * (Mpl_true[n] * MJtoMsun)^(0.3333) / ( rad_true[n] * RJtoAU );
        sini_onethird <- (1 - cosi[n]*cosi[n])^(0.16667);

	if ((x[n]/sini_onethird > xl) && (x[n]/sini_onethird < xu))
        {
		sum_log_x <- sum_log_x + log((gamma-ind[n]) / (xu^(gamma-ind[n]) - xl^(gamma-ind[n]))) + (gamma-1.-ind[n]) * log(x[n]/sini_onethird);
        }
        else
        {
	        sum_log_x <- -1/0.;
        }
    }

    for (n in 1:Ni)
    {
	xi[n] <- 0.462 * (peri_true[n]/365.242)^(0.66667) * (Mpli_true[n] * MJtoMsun)^(0.3333) / ( radi_true[n] * RJtoAU );

        if ((xi[n] > xl) && (xi[n] < xu))
        {
	        sum_log_x <- sum_log_x + log((gamma-indi[n]) / (xu^(gamma-indi[n]) - xl^(gamma-indi[n]))) + (gamma-1.-indi[n]) * log(xi[n]);
        }
    	else
        {
	        sum_log_x <- -1/0.;
        }
    }
        
    // specify priors in population papers
    gamma ~ uniform(-10.,10.);
    xl ~ uniform(0.,30.);
    xu ~ uniform(xl,30.);
    
    increment_log_prob( sum_log_x );
}
