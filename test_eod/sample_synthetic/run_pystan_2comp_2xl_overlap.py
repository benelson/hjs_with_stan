import numpy as np
import csv
import getData
import generate
import pystan

filename = "planets_synthetic.in"
modelfile = "powerlaw_2comp_2xl_overlap.stan"
nChains = 10
nIterations = 200000
nThin = 100
nJobs = -1

xi_data = [float(line.strip()) for line in open(filename, 'r')]

Stan_data = { 'Ni': len(xi_data), 'K': 2, 'xi': xi_data, 'x':[], 'N': 0}
Stan_ICs = generate.ICs_2comp_2xl_overlap([[0.25,0.75], [1.,-1.], [0.95, 2.], 21.], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_2xl_overlap(fit.extract())

with open("postSamp_synthetic_2comp_2xl_overlap.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
