data {
    int<lower=1> K; // number of mixture components
    int<lower=0> N; // number of planets w/o incl measurements
    int<lower=0> Ni; // number of planets w/ incl measurements
    real x[N]; // x values for N planets w/o incl measurements
    real xi[Ni]; // x values for Ni planets w/ incl measurements
}
parameters {
    simplex[K] theta; // mixing proportions
    real gamma[K]; // power-law indices
    real xl; // xlower
    real xt; // xtransition
    real xu; // xupper
    real cosi[N]; // cos inclinations of all planets
}
model {
    real ps[K];
    real xlower[K];
    real xupper[K];
    real integral[K];
    real sini_onethird;
    
    xlower[1] <- xl; // xl for first component
    xlower[2] <- xt; // xl for second component
    xupper[1] <- xt; // xu for first component
    xupper[2] <- xu; // xu for second component
    
    gamma ~ uniform(-20.,20.); // prior on gamma
    xt ~ uniform(0., 999.); // prior on xt
    xl ~ uniform(0.,xt); // prior on xl
    xu ~ uniform(xt, 999.); // prior on xu
    cosi ~ uniform(-1.,1.); // prior in cosi

    // for systems without inclination measurements
    for (n in 1:N) { //loop over planets
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
