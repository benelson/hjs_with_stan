data {
    int<lower=1> K; // number of mixture components
    int<lower=0> Ni; // number of planets w/ incl measurements
    
    real<lower=0.> peri[Ni];
    real<lower=0.> radi[Ni];
    real<lower=0.> Mpli[Ni];

    real<lower=0.> per_sigmai[Ni];
    real<lower=0.> rad_sigmai[Ni];
    real<lower=0.> Mpl_sigmai[Ni];
}
parameters {
    simplex[K] theta; // mixing proportions
    real gamma[K]; // power-law indices
    real xl[K]; // xlower
    real xu; // xupper

    real<lower=0.> peri_true[Ni];
    real<lower=-3., upper=1.> log_radi_true[Ni];
    real<lower=-3., upper=2.> log_Mpli_true[Ni];
}
transformed parameters{
    real radi_true[Ni];
    real Mpli_true[Ni];

    for (i in 1:Ni)
    {
        radi_true[i] <- 10. ^ log_radi_true[i];
        Mpli_true[i] <- 10. ^ log_Mpli_true[i];
    }
}
model {
    real ps[K];
    real xlower[K];
    real xupper[K];
    real integral[K];
    real sini_onethird;

    real RJtoAU;
    real MJtoMsun;

    real xi[Ni];

    RJtoAU <- 0.000477894503;
    MJtoMsun <- 0.000954265748;
    
    xlower[1] <- xl[1]; // xl for first component
    xlower[2] <- xl[2]; // xl for second component
    xupper[1] <- xu; // xu for first component
    xupper[2] <- xu; // xu for second component
    
    gamma ~ uniform(-10.,10.); // prior on gamma
    xl ~ uniform(0.,30.); // prior on xl
    xu ~ uniform(xl[2], 30.); // prior on xu

    peri ~ normal(peri_true, per_sigmai);
    radi ~ normal(radi_true, rad_sigmai);
    Mpli ~ normal(Mpli_true, Mpl_sigmai);

    for (n in 1:Ni) { //loop over planets
        xi[n] <- 0.462 * (peri_true[n]/365.242)^(0.66667) * (Mpli_true[n] * MJtoMsun)^(0.3333) / ( radi_true[n] * RJtoAU );
                
        for (k in 1:K) { //loop over mixture components
            if (n==1)
            {
                integral[k] <- (xupper[k]^(gamma[k]-1.) - xlower[k]^(gamma[k]-1.))/(gamma[k]-1.);
            }
            if ((xi[n] > xlower[k]) && (xi[n] < xupper[k]))
            {
                ps[k] <- theta[k] * (xi[n])^(gamma[k]-2.) / integral[k];
            }
            else
            {
                ps[k] <- 0.;
            }
        }     
        increment_log_prob( log(sum(ps)) );
    }
}
