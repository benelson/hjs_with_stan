import numpy as np
import csv
import generate
import getData
import pystan

def combine_two_dicts(x, y):
    z = {}
    for k in set(x).intersection(y):
            z[k] = np.array(list(x[k]) + list(y[k]))
    return z

filename_Kepler = "TableKeplerGiants.tex"
filename_rv = "planets_rvall.in"
filename_HAT = "planets_HAT.in"
filename_WASP = "planets_WASP.in"
modelfile = "powerlaw_2comp_overlap_at_unknown.stan"
nChains = 20
nIterations = 1000000
nThin = 1000
nJobs = -1

data_Kepler = getData.Santerne_to_dict(filename_Kepler)
data_rv = getData.exoplanetsORG_to_dict(filename_rv)
data_HAT = getData.exoplanetsORG_to_dict(filename_HAT)
data_WASP = getData.exoplanetsORG_to_dict(filename_WASP)

data_tmp = combine_two_dicts(data_Kepler,data_rv)
data_tmp = combine_two_dicts(data_tmp, data_HAT)
data_all = combine_two_dicts(data_tmp, data_WASP)

Stan_data = getData.create_Stan_input_mixture(data_all)
print(Stan_data)

Stan_ICs = generate.ICs_2comp_overlap_unknown([[0.25,0.75], [0.,-2.], 2.3, 25.], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_overlap_unknown(fit.extract())

with open("postSamp_all_2comp_overlap_at_unknown.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
