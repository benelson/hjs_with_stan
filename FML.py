import numpy as np

import compute
import create
import RVmodel

from scipy.stats import truncnorm
from math import erf

from collections import OrderedDict


class FML:
    samples = None
    obs = None
    nImportSamp = 0
    scale = 0.

    logAvg = 0.
    f_MCMC = 0.
    logFML = 0.


def lnlike_1comp(theta, obs):
    gamma, xl, xu = theta
    n = len(obs)
    norm = n * np.log(gamma / (xu**gamma - xl**gamma) )
    prod = (gamma-1.) * sum([np.log(x) for x in obs ])
    return norm + prod

def lnprior(theta):
    '''
    The log-prior function. Excluding prior in offset since we aren't comparing models with different offsets.
        
         Parameters
         ----------
         theta : 1-d list of Keplerian orbital parameters (p, K, e, w, M)
                 for every modeled planet
    '''
    gamma, xl, xu = theta

    lnp = 0.
    
    if -10. < gamma < 10. and 0. < xl < xu and xl < xu < 25.:
        lnp += -np.log(20.) - np.log(xu) - np.log(25. - xl)
        else:
            return -np.inf

    return lnp



def lnpost(theta, obs):
    '''
    The log-posterior function. Sum of log-likelihood and log-prior.
        
        Parameters
        ----------
        theta : 1-d list of Keplerian orbital parameters (p, K, e, w, M)
                for every modeled planet
        obs   : observations object with radial velocity measurement (.rv)
                and uncertainty (.err) attributes
    '''
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, obs)


class computeFML(FML):
    def __init__(self, samples, obs, nImportSamps=10000, scale=10.0):
        """
        Computes fully marginalized likelihood

        Parameters
        ----------
        samples      : posterior samples from emcee
        nPlanets     : number of planets
        nImportSamps : number of importance samples
        scale        : sets the scale for the truncation of the multivariate normal, measured in "sigmas"
        """

        self.samples = samples
        self.obs = obs
        self.nImportSamps = nImportSamps
        self.scale = scale

        param_keys, param_IS_keys = create.dict_keys(self.nPlanets, self.nOffsets, slope=self.slope)
        print(param_keys)
        print(param_IS_keys)
        postSamp, nPostSamples = create.posterior_samples_from_emcee(self.samples, param_keys)

        postSamp_pKhkl = compute.pKewM_to_importSamp_parameterization(postSamp, param_IS_keys, self.nPlanets)

        self.mediansG, self.covMatrixG, self.choleskyDecomp, self.logDetSigmaG = compute.matrix_info(postSamp_pKhkl)

        nParams = len(param_IS_keys)
        random_values = [ truncnorm.rvs(-self.scale, self.scale, size=nParams) for i in range(self.nImportSamps) ]

        samples = [ [] for i in range(self.nImportSamps) ]
        g_samples = [ [] for i in range(self.nImportSamps) ]
        loggs = [ 0. for i in range(self.nImportSamps) ]

        print("## Drawing importance samples...")

        for x in range(self.nImportSamps):
            dispersion = np.dot( self.choleskyDecomp, np.transpose(random_values[x]) )
            samples[x] = self.mediansG + dispersion
            g_samples[x] = list(samples[x])

            diff = np.subtract(samples[x],self.mediansG)

            logg = -0.5 * (nParams*np.log(2.*np.pi) + self.logDetSigmaG + \
                    np.dot( np.transpose(diff), \
                    np.linalg.solve(self.covMatrixG, np.subtract(samples[x],self.mediansG) ) ) ) - \
                    nParams*np.log(erf(self.scale/np.sqrt(2.)))
            loggs[x] = logg

        print("## Done drawing importance samples!")
        print("")
    
        g_samples_T = np.transpose(g_samples)
        importSamp_dict = OrderedDict()

        for i, item in enumerate(g_samples_T):
            importSamp_dict[param_IS_keys[i]] = item

        importSamp_pKhkl_dict = compute.importSamp_parameterization_to_pKewM(importSamp_dict, param_keys, self.nPlanets, self.pRatio)
        importSamp_pKewM = np.transpose([ vals for key, vals in importSamp_pKhkl_dict.items() ])

        print("## Evaluating lnpost at importance samples...")

        logPosteriors = np.array([ np.nan for i in range(self.nImportSamps) ])
        for i in range(nImportSamps):
            logPosteriors[i] = lnpost(importSamp_pKewM[i], obs, self.nPlanets, slope=self.slope)

        print("## Done evaluating lnpost!")
        print("")


        logSum = -(9.**99.)

        for i in range(self.nImportSamps):    
            diff = logPosteriors[i] - loggs[i]

            logSum = np.logaddexp(logSum, diff)
            if i%1000==0:
                print(str(i+1) + " " + str(logSum - np.log(i+1)))
    
        self.logAvg = logSum - np.log(self.nImportSamps)
        self.f_MCMC = 0.

        print("")
        print("## logAvg: " + str(self.logAvg))

        postSamp_wo_keys = []
        for key in postSamp_pKhkl:
            postSamp_wo_keys.append(postSamp_pKhkl[key])
    
        postSamp_wo_keys = np.transpose(np.array(postSamp_wo_keys))
        diff = postSamp_wo_keys-self.mediansG

        for j in range(nPostSamples):

            z = np.linalg.solve(self.choleskyDecomp, diff[j])

            if all([abs(k)<=scale for k in z]):
                self.f_MCMC += 1.
            else:
                self.f_MCMC += 0.
        
        self.f_MCMC = self.f_MCMC/nPostSamples
        self.logFML = self.logAvg - np.log(self.f_MCMC)

        print("## f_MCMC: " + str(self.f_MCMC))
        print("## logFML: " + str(self.logFML))

        print("## FML computed!")
        print("## Done!")
