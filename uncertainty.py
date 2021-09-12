#uncertainty
#contains functions for uncertainty from lab 4

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from solve_ode import solve_ode

text_size = 16
np.random.seed(69) #Nice

def grid_search(b_1, b_2, b_3, tau, p_0, m_0, alpha):

    #get data
    t_calibrate, nitrate_calibrate = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)

    # number of values considered for each parameter within a given interval
    N = 10

    # vectors of parameter values
    b_1 = np.linspace(b_1*0.95,b_1/0.95, N)
    b_2 = np.linspace(b_2*0.95,b_2/0.95, N)

    # grid of parameter values: returns every possible combination of parameters in b_1 and b_2
    B_1, B_2 = np.meshgrid(b_1, b_2, indexing='ij')

    # empty 2D matrix for objective function
    S = np.zeros(B_1.shape)

    # error variance - 4 mg/L
    v = 4.

    #get t values for interpolation
    t_numeric, n_numeric, _ = solve_ode(b_1=b_1[1], b_2=b_2[1], b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast= True)
    
    #interpolate nitrate_calibrate values to align with t_numeric
    nitrate_cal_int = np.interp(t_numeric, t_calibrate, nitrate_calibrate)

    # grid search algorithm
    for i in range(len(b_1)):
        for j in range(len(b_2)):

    #compute the sum of squares objective function at each value 
            t_numeric, n_numeric, _ = solve_ode(b_1=b_1[i], b_2=b_2[i], b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast = True)
            S[i,j] = np.sum((nitrate_cal_int-n_numeric)**2)/v

    #compute the posterior
    P = np.exp(-S/2.)
    Pint = np.sum(P)*(b_1[1]-b_1[0])*(b_2[1]-b_2[0])
    P = P/Pint

    return b_1,b_2,P

def construct_samples(b1,b2,P,N_samples,b_1, b_2, b_3, tau, p_0, m_0, alpha):
    ''' This function constructs samples from a multivariate normal distribution
        fitted to the data.
        Parameters:
        -----------
        b1 : array-like
            Vector of 'a' parameter values.
        b2 : array-like
            Vector of 'b' parameter values.
        P : array-like
            Posterior probability distribution.
        N_samples : int
            Number of samples to take.
        Returns:
        --------
        samples : array-like
            parameter samples from the multivariate normal
    '''

    # compute properties (fitting) of multivariate normal distribution
    # mean = a vector of parameter means
    # covariance = a matrix of parameter variances and correlations
    A, B = np.meshgrid(b1,b2,indexing='ij')
    mean, covariance = fit_mvn([A,B], P)

    #create samples using numpy function multivariate_normal
    samples = np.random.multivariate_normal(mean, covariance, N_samples)

    return samples

def model_ensemble(samples,b_1, b_2, b_3, tau, p_0, m_0, alpha):
    ''' Runs the model for given parameter samples and plots the results.
        Parameters:
        -----------
        samples : array-like
            parameter samples from the multivariate normal
    '''

    # get the data
    tp,po = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)	

    #multipliers for cattle population used as input to solve_ode
    cattle_multipliers = [0, 0.5, 1, 2]

    #values for MAR pressure
    mar_values = [0, 0.1, 0.3]

    #nested for loop to plot all samples for all scenarios
    for k in mar_values:

        #create figure, axis, for each mar value plot
        fig, ax = plt.subplots()
        ax.set_xlabel('time (yrs)')
        ax.set_ylabel('nitrate concentration (mg/L)')
        ax.set_title(f'Forecast of Concentration of Nitrate in Southland Aquifer with MAR at a magnitude of {k} MPa')

        #filter through cattle populations
        for m in cattle_multipliers:
            
            #plot individual samples
            for a,b in samples:
                t,pm,_ = solve_ode(a,b, b_3, tau, p_0, m_0, alpha, forecast = True, P_mar=k, multiplier=m)
                ax.plot(t,pm,'k-', lw=0.25,alpha=0.4)
            
            #plot best fit curves
            t_array, forecast, _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast=True, multiplier=m, P_mar=k)
            ax.plot(t_array, forecast, label=f"{m * 100}% of current cow population", lw = 2, alpha = 1)
        
        #plot data points over top
        ax.scatter(tp, po,c ='red', label = 'Nitrate Concentration Data')
        ax.plot([],[],'k-', lw=0.5,alpha=0.4, label='model ensemble')
        
        #generate legend. Ur a legend if ur reading this ;) 
        ax.legend()

        #used to save figure
        #plt.savefig(f"forecast_{j}.jpg")

        #show plot
        plt.show()


def fit_mvn(parspace, dist):
    """Finds the parameters of a multivariate normal distribution that best fits the data
    Parameters:
    -----------
        parspace : array-like
            list of meshgrid arrays spanning parameter space
        dist : array-like 
            PDF over parameter space
    Returns:
    --------
        mean : array-like
            distribution mean
        cov : array-like
            covariance matrix		
    """

    # dimensionality of parameter space
    N = len(parspace)

    # flatten arrays
    parspace = [p.flatten() for p in parspace]
    dist = dist.flatten()

    # compute means
    mean = [np.sum(dist*par)/np.sum(dist) for par in parspace]

    # compute covariance matrix
        # empty matrix
    cov = np.zeros((N,N))
        # loop over rows
    for i in range(0,N):
            # loop over upper triangle
        for j in range(i,N):
                # compute covariance
            cov[i,j] = np.sum(dist*(parspace[i] - mean[i])*(parspace[j] - mean[j]))/np.sum(dist)
                # assign to lower triangle
            if i != j: cov[j,i] = cov[i,j]
            
    return np.array(mean), np.array(cov)




