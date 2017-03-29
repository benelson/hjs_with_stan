import numpy as np
import csv
import generate
import getData
import pystan
import sys

arg = sys.argv[1]
filename = "planets_HAT.csv"
modelfile = "powerlaw_2comp_2xl_overlap_alltransit.stan"
nChains = 1
nIterations = 400000
nThin = 200
nJobs = 1

data = getData.NEA_to_dict(filename)
Stan_data = getData.create_Stan_input_mixture_alltransit(data)

Stan_ICs = generate.ICs_2comp_2xl_overlap_alltransit([[0.35,0.65], [0.5,-1.5], [1.35, 2.3], 18.], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_2xl_overlap(fit.extract())

with open("postSamp_HAT_2comp_2xl_overlap." + arg + ".txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))

