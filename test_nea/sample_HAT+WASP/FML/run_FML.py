import FML
import compute
import create
import observations

obs_file = system+"_data"

post_samp_file = "../postSamp_HAT+WASP_2comp_2xl_overlap.txt"
obs = observations.get_obs(obs_file)
post_samps = create.post_samp_from_file(post_samp_file)

FML_data = FML.computeFML(post_samps, obs=obs, nPlanets=npl, nOffsets=noff, nImportSamps=300000, scale=10.0, pRatio=tmp, slope=slope)
