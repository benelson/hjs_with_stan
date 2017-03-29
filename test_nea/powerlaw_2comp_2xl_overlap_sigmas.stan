data {
    int<lower=1> K; // number of mixture components
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
}
parameters {
    simplex[K] theta; // mixing proportions
    real gamma[K]; // power-law indices
    real xl[K]; // xlower
    real xu; // xupper
    real cosi[N]; // cos inclinations of all planets

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
    real ps[K];
    real xlower[K];
    real xupper[K];
    real integral[K];
    real sini_onethird;

    real RJtoAU;
    real MJtoMsun;

    real x[N];
    real xi[Ni];

    RJtoAU <- 0.000477894503;
    MJtoMsun <- 0.000954265748;
    
    xlower[1] <- xl[1]; // xl for first component
    xlower[2] <- xl[2]; // xl for second component
    xupper[1] <- xu; // xu for first component
    xupper[2] <- xu; // xu for second component
    
    gamma ~ uniform(-20.,20.); // prior on gamma
    xl ~ uniform(0.,25.); // prior on xl
    xu ~ uniform(xl[2], 25.); // prior on xu
    cosi ~ uniform(-1.,1.); // prior in cosi

//    per_true ~ uniform(0.1,100.);
//    peri_true ~ uniform(0.1,100.);

    per ~ normal(per_true, per_sigma);
    rad ~ normal(rad_true, rad_sigma);
    Mpl ~ normal(Mpl_true, Mpl_sigma);

    peri ~ normal(peri_true, per_sigmai);
    radi ~ normal(radi_true, rad_sigmai);
    Mpli ~ normal(Mpli_true, Mpl_sigmai);


    // for systems without inclination measurements
    for (n in 1:N) { //loop over planets
        x[n] <- 0.462 * (per_true[n]/365.242)^(0.66667) * (Mpl_true[n] * MJtoMsun)^(0.3333) / ( rad_true[n] * RJtoAU );
        sini_onethird <- (1 - cosi[n]*cosi[n])^(0.16667); // converting cosi to sini^(1/3)
        
        for (k in 1:K) { //loop over mixture components
            if (n==1)
            {
                integral[k] <- (xupper[k]^gamma[k] - xlower[k]^gamma[k])/gamma[k];
            }
            if ((x[n]/sini_onethird > xlower[k]) && (x[n]/sini_onethird < xupper[k]))
            {
                ps[k] <- theta[k] * (x[n]/sini_onethird)^(gamma[k]-1.) / integral[k];
            }
            else
            {
                ps[k] <- 0.;
            }
        }     
        increment_log_prob( log(sum(ps)) );
    }
    
    for (n in 1:Ni) { //loop over planets
        xi[n] <- 0.462 * (peri_true[n]/365.242)^(0.66667) * (Mpli_true[n] * MJtoMsun)^(0.3333) / ( radi_true[n] * RJtoAU );
                
        for (k in 1:K) { //loop over mixture components
            if (n==1)
            {
                integral[k] <- (xupper[k]^gamma[k] - xlower[k]^gamma[k])/gamma[k];
            }
            if ((xi[n] > xlower[k]) && (xi[n] < xupper[k]))
            {
                ps[k] <- theta[k] * (xi[n])^(gamma[k]-1.) / integral[k];
            }
            else
            {
                ps[k] <- 0.;
            }
        }     
        increment_log_prob( log(sum(ps)) );
    }
}
