import numpy as np

def ICs_1comp(IC, xlist, nChains):
    gamma = IC[0]+np.random.randn(nChains)*0.05
    xl = IC[1]+np.random.randn(nChains)*0.01
    xu = IC[2]+np.random.randn(nChains)*0.01

    ICs = []
    for i in range(nChains):
        ICs.append({ "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs

def ICs_2comp_1xl(IC, xlist, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xl = [ IC[2]*0.95 for i in range(nChains) ]
    xu = [ IC[3]*1.05 for i in range(nChains) ]

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs

def ICs_2comp_2xl(IC, xlist, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xl = [ IC[2]*0.95 for i in range(nChains) ]
    xt = [ IC[3] for i in range(nChains) ]
    xu = [ IC[4]*1.05 for i in range(nChains) ]
    
    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xt": xt[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs


def ICs_2comp_0xl_overlap(IC, xlist, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xu = [ IC[2] for i in range(nChains) ]

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs


def ICs_2comp_1xl_overlap(IC, xlist, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xl = [ IC[2] + np.random.normal(0.,0.1) for i in range(nChains) ]
    xu = [ IC[3] for i in range(nChains) ]

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs


def ICs_2comp_2xl_overlap(IC, xlist, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xl = [ IC[2] + np.random.randn(2)*0.05 for i in range(nChains) ]
    xu = [ IC[3] for i in range(nChains) ]

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs



def ICs_1comp_sigmas(IC, data, nChains):
    gamma = IC[0]+np.random.randn(nChains)*0.05
    xl = IC[1]+np.random.randn(nChains)*0.01
    xu = IC[2]+np.random.randn(nChains)*0.01

    per = data["per"]
    log_rad = list(np.log10(data["rad"]))
    log_Mpl = list(np.log10(data["Mpl"]))
#    Mstar = data["Mstar"]

    peri = data["peri"]
    log_radi = list(np.log10(data["radi"]))
    log_Mpli = list(np.log10(data["Mpli"]))
#    Mstari = data["Mstari"]

    N = len(per)

    ICs = []
    for i in range(nChains):
        ICs.append({ "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(N)],\
                     "per_true": per, "log_rad_true": log_rad, "log_Mpl_true": log_Mpl ,\
                     "peri_true": peri, "log_radi_true": log_radi, "log_Mpli_true": log_Mpli })
    return ICs


def ICs_2comp_2xl_sigmas(IC, data, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xl = [ IC[2]*0.95 for i in range(nChains) ]
    xt = [ IC[3] for i in range(nChains) ]
    xu = [ IC[4]*1.05 for i in range(nChains) ]

    per = data["per"]
    log_rad = list(np.log10(data["rad"]))
    log_Mpl = list(np.log10(data["Mpl"]))

    peri = data["peri"]
    log_radi = list(np.log10(data["radi"]))
    log_Mpli = list(np.log10(data["Mpli"]))

    N = len(per)

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xt": xt[i], "xu": xu[i], "cosi": [0. for i in range(N)],\
                         "per_true": per, "log_rad_true": log_rad, "log_Mpl_true": log_Mpl ,\
                         "peri_true": peri, "log_radi_true": log_radi, "log_Mpli_true": log_Mpli })
    return ICs


def ICs_2comp_2xl_overlap_sigmas(IC, data, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.04 for i in range(nChains) ]
    xl = [ IC[2] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xu = [ IC[3]*1.05 for i in range(nChains) ]

    per = data["per"]
    log_rad = list(np.log10(data["rad"]))
    log_Mpl = list(np.log10(data["Mpl"]))

    peri = data["peri"]
    log_radi = list(np.log10(data["radi"]))
    log_Mpli = list(np.log10(data["Mpli"]))

    N = len(per)

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(N)],\
                         "per_true": per, "log_rad_true": log_rad, "log_Mpl_true": log_Mpl ,\
                         "peri_true": peri, "log_radi_true": log_radi, "log_Mpli_true": log_Mpli })
    return ICs
