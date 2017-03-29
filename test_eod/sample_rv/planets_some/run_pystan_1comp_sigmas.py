import csv
import generate
import getData
import pystan

filename = "planets_rvsome.in"
modelfile = "powerlaw_1comp_sigmas.stan"
nChains = 10
nIterations = 200000
nThin = 100
nJobs = -1

data = getData.exoplanetsORG_to_dict(filename)
Stan_data = getData.create_Stan_input_sigmas(data, transit=0)
Stan_ICs = generate.ICs_1comp_sigmas([-1.,1.5,30.], Stan_data, nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_1comp(fit.extract())

print(postSamp)

with open("postSamp_rvsome_1comp_sigmas.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
