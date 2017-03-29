import numpy as np
import csv
import generate
import getData
import pystan
import sys

def combine_two_dicts(x, y):
    z = {}
    for k in set(x).intersection(y):
            z[k] = np.array(list(x[k]) + list(y[k]))
    return z

filename_HAT = "planets_HAT.in"
filename_WASP = "planets_WASP.in"
modelfile = "powerlaw_1comp_alltransit.stan"
nChains = 1
nIterations = 200000
nThin = 100
nJobs = 1

arg = sys.argv[1]

data_HAT = getData.exoplanetsORG_to_dict(filename_HAT)
data_WASP = getData.exoplanetsORG_to_dict(filename_WASP)
data_all = combine_two_dicts(data_HAT,data_WASP)

Stan_data = getData.create_Stan_input_alltransit(data_all)
print(Stan_data)
Stan_ICs = generate.ICs_1comp_alltransit([-1., 1.05, 13.], nChains)
print(Stan_ICs)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_1comp(fit.extract())

with open("postSamp_HAT+WASP_1comp." + arg + ".txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
