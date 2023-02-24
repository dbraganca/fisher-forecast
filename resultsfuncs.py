import numpy as np
import matplotlib.pyplot as plt
from getdist import plots, MCSamples
from scipy.io import mmread
from tabulate import tabulate

def multiindex(arr1, arr2):
    """
    This function takes in two lists, `arr1` and `arr2`, and returns a new list containing the indices in `arr1` 
    that match the elements in `arr2`. 
    
    :param arr1: list, a list of elements to search through
    :param arr2: list, a list of elements to search for
    :return: list, a list of indices in `arr1` that match the elements in `arr2`
    """
    res = []
    for i in range(len(arr2)):
        new = list(np.where(np.array(arr1) == arr2[i])[0])   
        res.extend(new)
    return list(np.array(res).flatten())

def mutransform(cov,mupos,backm):
    """
    Fisher matrix transformation from log neutrino mass to neutrino mass
    """
    cov[mupos,:]=cov[mupos,:]*backm
    cov[:,mupos]=cov[:,mupos]*backm
    return cov

# build parameters
cosmoindex = ["lnAs", "ns", "h", "Om", "ob", "mn", "Ok"]

fnlindex = ["fNLloc", "fNLeq", "fNLorth"]

biasindex = ["b1", "c2", "b3", "b4", "c4", "b6",
             "b7", "b8", "b9", "b10", "b11", "Bc1", 
             "Bc2", "Bc3", "Bc4", "Be1", "Be2", "ce2", 
             "Bc5", "Bc6", "Bc7", "Bc8", "Bc9", "Bc10", 
             "Bc11", "Bc12", "Bc13", "Bc14", "Bd1", 
             "Bd2", "Bd3", "Be3", "Be4", "Be5",
             "Be6", "Be7", "Be8", "Be9", "Be10","Be11", "Be12", "Tst"]

paramindex = cosmoindex + fnlindex + biasindex


# paramindex = ["lnAs", "ns", "h", "Om", "ob", "mn",
# "fNLloc", "fNLeq", "fNLorth", "b1", "c2", "b3", "b4", "c4", "b6",
# "b7", "b8", "b9", "b10", "b11", "Bc1", "Bc2", "Bc3", "Bc4", "Be1",
# "Be2", "ce2", "Bc5", "Bc6", "Bc7", "Bc8", "Bc9", "Bc10", "Bc11",
# "Bc12", "Bc13", "Bc14", "Bd1", "Bd2", "Bd3", "Be3", "Be4", "Be5",
# "Be6", "Be7", "Be8", "Be9", "Be10","Be11", "Be12", "Tst"]

# bias list without Be1 and Bd1
biasindex2 = ["b1", "c2", "b3", "b4", "c4", "b6",
              "b7", "b8", "b9", "b10", "b11", "Bc1", "Bc2", "Bc3", "Bc4",
              "Be2", "ce2", "Bc5", "Bc6", "Bc7", "Bc8", "Bc9", "Bc10", "Bc11",
              "Bc12", "Bc13", "Bc14", "Bd2", "Bd3", "Be3", "Be4", "Be5",
              "Be6", "Be7", "Be8", "Be9", "Be10","Be11", "Be12", "Tst"]

def parampos(arr): 
    return multiindex(paramindex,arr)


paramindex4sky = cosmoindex + fnlindex + list(
                    np.transpose(np.array(
                        [np.array(biasindex),np.array(biasindex),np.array(biasindex),np.array(biasindex)]
                        )).flatten())

num_param_1sky = len(paramindex)
num_param_4sky = len(paramindex4sky)


def parampos4sky(arr): 
    return multiindex(paramindex4sky,arr)


Paramlabels = ["ln(10^{10}A_s)", "n_s", "h", "\Omega_m", "\omega_b", "\sum_im_{\\nu_i}", "\Omega_k",
               "f_{NL}^{loc}", "f_{NL}^{eq}", "f_{NL}^{orth}", 
               "b_1", "c_2", "b_3", "b_4", "c_4", "b_6", "b_7", "b_8", "b_9",
               "b_{10}", "b_{11}", "Bc_1", "Bc_2", "Bc_3", "Bc_4","Be_1", 
               "Be_2", "ce_2", "Bc_5", "Bc_6", "Bc_7", "Bc_8", "Bc_9", "Bc_{10}",
               "Bc_{11}", "Bc_{12}", "Bc_{13}", "Bc_{14}", "Bd_1", "Bd_2", "Bd_3", "Be_3", "Be4",
               "Be_5", "Be_6", "Be_7", "Be_8", "Be_9", "Be_{10}", "Be_{11}", "Be_{12}", "Tst"]

background = [3.044, 0.965, 0.673, 0.322106, 0.02237, 0.3,0,0,0,1.9391, 
            1.1447, -0.374, 0.1276, -0.2895, -0.3487, 0.2234, -0.2972, 0.0148, 
            0.0428, 0.0356, 5.4605, -1.5439, 1.3081, -0.4785, 1.6922, 0.9091, 
            0.5548, 0.1054, 0.8713, -0.4555, 0.4442, -0.4152, -0.65, -0.0881, 
            -0.3731, -0.1642, -0.1958, 5.4, -0.721, -0.4266, 0.0737, -0.1381, 
            6.0673, -0.0932, -0.9702, 0.2617, 0.2638, -0.1526, 0.427, -0.4322, 
            9.6]

def geterror(cov):
    """
    Compute the standard errors from a covariance matrix.

    Parameters:
    cov (numpy.ndarray): A square covariance matrix.

    Returns:
    numpy.ndarray: A 1D array of standard errors.
    """
    errors = []
    for i in range(len(cov)):
        errors.append(np.sqrt(cov[i,i]))
    return errors

def around2(arr,n):
    """
    Round a numpy array to a given number of decimal places.

    Parameters:
    arr (numpy.ndarray): A numpy array to be rounded.
    n (int): The number of decimal places to round to.

    Returns:
    numpy.ndarray: A numpy array with the values rounded to n decimal places.
    """
    rounded = np.zeros(arr.shape)
    for i in range(len(arr)):
        for j in range(len(arr[0])):
            rounded[i,j]=np.format_float_positional(arr[i,j], precision=n, unique=False, fractional=False, trim='k')
    return rounded


def gettable(errors,cosmopar,labs):
    """
    Create a table of cosmological parameter errors.

    Parameters:
    errors (list): A list of standard errors for each parameter.
    cosmopar (list): A list of cosmological parameter names.
    labs (list): A list of labels for the parameter names.

    Returns:
    list: A nested list representing the table.
    """
    errtab = errors
    for i in range(len(labs)):
        errtab[i].insert(0,labs[i])
    errtab = np.array(errtab).T.tolist()
    errtab[0].insert(0,r"\sigma")
    for i in range(len(cosmopar)):
        errtab[i+1].insert(0,cosmopar[i])
    return np.array(errtab).T.tolist()

def combinesamps(samplist):
    """
    Combine a list of MCMC samples into a single sample.

    Parameters:
    samplist (list): A list of MCMC samples.

    Returns:
    numpy.ndarray: A combined array of samples.
    """
    samps = []
    for j in range(len(samplist)):
        for i in range(len(samplist[j])):
            samps.append(samplist[j][i])
    return samps

def getPlot(base,cosmovar,plotvar,labs,save=0):
    """
    Plot a triangle plot of cosmological parameter constraints.

    Parameters:
    base (list): A list of covariance matrices.
    cosmovar (list): A list of cosmological parameter names.
    plotvar (list): A list of parameter names to plot.
    labs (list): A list of labels for the parameter names.
    save (int): Whether to save the plot and parameter table to files (default: 0).

    Returns:
    None
    """
    samps = []
    errors = []
    cosmopar = list(np.array(Paramlabels)[parampos(plotvar)])
    for i in range(len(base)):
        if(len(base[i])==num_param_4sky):
            freeparams = parampos4sky(cosmovar)
            newparams =list(np.array(paramindex4sky)[freeparams])
            Plotparams = multiindex(newparams,plotvar)    
        elif(len(base[i])==num_param_1sky):
            freeparams = parampos(cosmovar)
            newparams =list(np.array(paramindex)[freeparams])
            Plotparams = multiindex(newparams,plotvar)
        else:
            print('wrong dim')
            
        Fisher = base[i][freeparams,:][:,freeparams]
        if 'mn' in cosmovar:
            mupos = list(np.where(np.array(cosmovar) == 'mn')[0])[0]  
            muback= background[5]
            fishcov = mutransform(np.linalg.inv(Fisher),mupos,muback)[Plotparams,:][:,Plotparams]
        else:
            fishcov = np.linalg.inv(Fisher)[Plotparams,:][:,Plotparams]
        sampl = np.random.multivariate_normal(list(np.array(background)[parampos(plotvar)]), fishcov, size=int(1e5))
        samps.append(MCSamples(samples=sampl, names=cosmopar, labels = cosmopar))
        errors.append(geterror(fishcov))
        
    g = plots.get_subplot_plotter()
    g.triangle_plot(samps,cosmopar, filled=True, legend_labels = labs,title_limit=1)
    if(save):
        plt.savefig(save+'.pdf')
        with open(save+'.txt', 'w') as f:
            f.write(tabulate(gettable(around2(np.array(errors),2).tolist(),cosmopar,labs),tablefmt="latex")) 

def getSamps(base,cosmovar,plotvar,labs):

    samps = []
    cosmopar = list(np.array(Paramlabels)[parampos(plotvar)])

    for i in range(len(base)):
        if(len(base[i])==num_param_4sky):
            freeparams = parampos4sky(cosmovar)
            newparams =list(np.array(paramindex4sky)[freeparams])
            Plotparams = multiindex(newparams,plotvar)    
        elif(len(base[i])==num_param_1sky):
            freeparams = parampos(cosmovar)
            newparams =list(np.array(paramindex)[freeparams])
            Plotparams = multiindex(newparams,plotvar)
        else:
            print('wrong dim')
            
        Fisher = base[i][freeparams,:][:,freeparams]
        if 'mn' in cosmovar:
            mupos = list(np.where(np.array(cosmovar) == 'mn')[0])[0]  
            muback= background[5]
            fishcov = mutransform(np.linalg.inv(Fisher),mupos,muback)[Plotparams,:][:,Plotparams]
        else:
            fishcov = np.linalg.inv(Fisher)[Plotparams,:][:,Plotparams]
        sampl = np.random.multivariate_normal(list(np.array(background)[parampos(plotvar)]), fishcov, size=int(1e5))
        samps.append(MCSamples(samples=sampl, names=cosmopar, labels = cosmopar))
    return samps ,cosmopar,labs


def getPlotfromSamps(samps,cosmopar,labs,save=0):
    g = plots.get_subplot_plotter()
    g.triangle_plot(samps,cosmopar, filled=True, legend_labels = labs,title_limit=1)
    if(save):
        plt.savefig(save+'.pdf')

def getError(base,cosmovar,plotvar):
    """
    Calculates errors for each element in 'base' by computing the inverse of the Fisher matrix,
    selecting the relevant parameters defined by 'plotvar', and computing the error on those parameters.

    Args:
    - base (list): A list of numpy arrays, where each array is a Fisher matrix. The Fisher matrix
                   corresponds to a specific set of cosmological parameters.
    - cosmovar (list): A list of strings, where each string represents a cosmological parameter. This list
                       defines the order of the parameters in the Fisher matrix.
    - plotvar (list): A list of integers, where each integer represents a specific combination of
                      cosmological parameters. This list defines which parameters to plot in the output.

    Returns:
    - errors (list): A list of errors for each element in 'base'. The error is calculated by computing
                     the inverse of the Fisher matrix, selecting the relevant parameters defined by 'plotvar',
                     and computing the error on those parameters.
                     
    Note:
    This function makes use of the numpy and getDist packages.
    """

    samps = []
    errors = []
    cosmopar = list(np.array(Paramlabels)[parampos(plotvar)])
    for i in range(len(base)):
        if(len(base[i])==num_param_4sky):
            freeparams = parampos4sky(cosmovar)
            newparams =list(np.array(paramindex4sky)[freeparams])
            Plotparams = multiindex(newparams,plotvar)    
        elif(len(base[i])==num_param_1sky):
            freeparams = parampos(cosmovar)
            newparams =list(np.array(paramindex)[freeparams])
            Plotparams = multiindex(newparams,plotvar)
        else:
            print('wrong dim')
            
        Fisher = base[i][freeparams,:][:,freeparams]
        if 'mn' in cosmovar:
            mupos = list(np.where(np.array(cosmovar) == 'mn')[0])[0]  
            muback= background[5]
            fishcov = mutransform(np.linalg.inv(Fisher),mupos,muback)[Plotparams,:][:,Plotparams]
        else:
            fishcov = np.linalg.inv(Fisher)[Plotparams,:][:,Plotparams]
        sampl = np.random.multivariate_normal(list(np.array(background)[parampos(plotvar)]), fishcov, size=int(1e5))
        samps.append(MCSamples(samples=sampl, names=cosmopar, labels = cosmopar))
        errors.append(geterror(fishcov))
    return errors


def fnlerrors(base,cosmovar,plotvar,save=0):
    errors = []
    for i in range(len(base)):
        if(len(base[i])==num_param_4sky):
            freeparams = parampos4sky(cosmovar)
            newparams =list(np.array(paramindex4sky)[freeparams])
            Plotparams = multiindex(newparams,plotvar)    
        elif(len(base[i])==num_param_1sky):
            freeparams = parampos(cosmovar)
            newparams =list(np.array(paramindex)[freeparams])
            Plotparams = multiindex(newparams,plotvar)
        else:
            print('wrong dim')
        Fisher = base[i][freeparams,:][:,freeparams]
        if 'mn' in cosmovar:
            mupos = list(np.where(np.array(cosmovar) == 'mn')[0])[0]  
            muback= background[5]
            fishcov = mutransform(np.linalg.inv(Fisher),mupos,muback)[Plotparams,:][:,Plotparams]
        else:
            fishcov = np.linalg.inv(Fisher)[Plotparams,:][:,Plotparams]
        errors.append(np.sqrt(fishcov[0,0]))
        
    if (save):
        np.savetxt(save+'.txt',np.array(errors))
    print(errors)








