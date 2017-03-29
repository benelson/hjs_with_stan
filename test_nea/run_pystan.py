import csv
import generate
import getData
import pystan

dataset = "rv"

filename = "sample_"+dataset+"/planets.in"
modelfile = "powerlaw_1comp.stan"
nChains = 10
nIterations = 100000
nThin = 50
nJobs = 1

data = getData.exoplanetsORG_to_dict(filename)
Stan_data = getData.create_Stan_input(data)
Stan_ICs = generate.ICs_1comp([-1.,1.9,25.], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_1comp(fit.extract())

print(postSamp)

with open("postSamp"+dataset+".txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
