import csv
import generate
import getData
import pystan

filename = "planets_rvall.in"
modelfile = "powerlaw_1comp.stan"
nChains = 10
nIterations = 300000
nThin = 150
nJobs = -1

data = getData.exoplanetsORG_to_dict(filename)
Stan_data = getData.create_Stan_input(data, transit=0)
Stan_ICs = generate.ICs_1comp([-1.,1.5,25.], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_1comp(fit.extract())

with open("postSamp_rvall_1comp.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
