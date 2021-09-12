#uncertainty
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from solve_ode import solve_ode

    # 1. DEFINE parameter ranges for the grid search
    # 2. COMPUTE the sum-of-squares objective function for each parameter combination
    # 3. COMPUTE the posterior probability distribution
    # 4. ANSWER the questions in the lab document

    # 1. define parameter ranges for the grid search
text_size = 16
np.random.seed(69) #Nice

def plot_posterior2D(a, b, P):	
    """Plot posterior distribution

    Args:
        a (numpy array): a distribution vector
        b (numpy array): b distribution vector
        P (numpy array): posterior matrix
    """

    # grid of parameter values: returns every possible combination of parameters in a and b
    A, B = np.meshgrid(a, b)

    # plotting
    fig = plt.figure(figsize=[10., 7.])				# open figure
    ax1 = fig.add_subplot(111, projection='3d')		# create 3D axes
    ax1.plot_surface(A, B, P, rstride=1, cstride=1,cmap=cm.coolwarm, lw = 0.5,edgecolor='k')	# show surface

    # plotting upkeep
    ax1.set_xlabel('a', fontsize = text_size)
    ax1.set_ylabel('b', fontsize = text_size)
    ax1.set_zlabel('P', fontsize = text_size)
    ax1.set_xlim([a[0], a[-1]])
    ax1.set_ylim([b[0], b[-1]])
    ax1.set_zlim(0., )
    ax1.view_init(40, 100.)

    # save and show
    plt.show()

def grid_search(b_1, b_2, b_3, tau, p_0, m_0, alpha):

    #get data
    t_calibrate, nitrate_calibrate = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)

    # number of values considered for each parameter within a given interval
    N = 5

    # vectors of parameter values
    b_1 = np.linspace(b_1*0.95,b_1/0.95, N)
    b_2 = np.linspace(b_2*0.95,b_2/0.95, N)

    # grid of parameter values: returns every possible combination of parameters in b_1 and b_2
    B_1, B_2 = np.meshgrid(b_1, b_2, indexing='ij')

    # empty 2D matrix for objective function
    S = np.zeros(B_1.shape)

    # error variance - 2 bar
    v = 1.

    #get t values for interpolation 
    t_numeric, n_numeric, _ = solve_ode(b_1=b_1[1], b_2=b_2[1], b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast= True)
    nitrate_cal_int = np.interp(t_numeric, t_calibrate, nitrate_calibrate)
    # grid search algorithm
    for i in range(len(b_1)):
        for j in range(len(b_2)):
    # 3. compute the sum of squares objective function at each value 
    #pm =
    #S[i,j] =
            #pm = ode_solve(tp,a[i],b[j])
            t_numeric, n_numeric, _ = solve_ode(b_1=b_1[i], b_2=b_2[i], b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast = True)
            S[i,j] = np.sum((nitrate_cal_int-n_numeric)**2)/v

    # 4. compute the posterior
    #P=
    P = np.exp(-S/2.)
    Pint = np.sum(P)*(b_1[1]-b_1[0])*(b_2[1]-b_2[0])
    P = P/Pint
    plot_posterior2D(b_1,b_2,P)

    return b_1,b_2,P

def construct_samples(a,b,P,N_samples,b_1, b_2, b_3, tau, p_0, m_0, alpha):
    ''' This function constructs samples from a multivariate normal distribution
        fitted to the data.

        Parameters:
        -----------
        a : array-like
            Vector of 'a' parameter values.
        b : array-like
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
    # **to do**
    # 1. FIGURE OUT how to use the multivariate normal functionality in numpy
    #    to generate parameter samples
    # 2. ANSWER the questions in the lab document

    # compute properties (fitting) of multivariate normal distribution
    # mean = a vector of parameter means
    # covariance = a matrix of parameter variances and correlations
    A, B = np.meshgrid(a,b,indexing='ij')
    mean, covariance = fit_mvn([A,B], P)

    #create samples using numpy function multivariate_normal (Google it)
    #samples=
    samples = np.random.multivariate_normal(mean, covariance, N_samples)

    # plot samples and predictions
    plot_samples2D(a, b, P=P, samples=samples, b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha)

    return samples

def model_ensemble(samples,b_1, b_2, b_3, tau, p_0, m_0, alpha):
    ''' Runs the model for given parameter samples and plots the results.

        Parameters:
        -----------
        samples : array-like
            parameter samples from the multivariate normal
    '''
    # **to do**
    # Run your parameter samples through the model and plot the predictions.

    # 1. choose a time vector to evaluate your model between 1953 and 2012 
    # t =


    # 2. create a figure and axes (see TASK 1)
    #f,ax =
    f,ax = plt.subplots(1,1)

    # 3. for each sample, solve and plot the model  (see TASK 1)
    

    

    # get the data
    tp,po = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)	
    ax.axvline(1980, color='b', linestyle=':', label='calibration/forecast')

    # 4. plot Wairakei data as error bars
    # *hint* see TASK 1 for appropriate plotting commands
    v = 1.
    ax.errorbar(tp,po,yerr=v,fmt='ro', label='data')
    ax.set_xlabel('time')
    ax.set_ylabel('pressure')
    ax.legend()
    plt.show()


    cattle_multipliers = [0, 0.5, 1, 2]

    mar_values = [0, 0.1, 0.3]
    for k in mar_values:    
        fig, ax = plt.subplots()
        ax.set_xlabel('time (yrs)')
        ax.set_ylabel('nitrate concentration (mg/L)')
        ax.set_title(f'Forecast of Concentration of Nitrate in Southland Aquifer with MAR at a magnitude of {k} MPa')

        for m in cattle_multipliers:
            
            for a,b in samples:
                t,pm,_ = solve_ode(a,b, b_3, tau, p_0, m_0, alpha, forecast = True, P_mar=k, multiplier=m)
                ax.plot(t,pm,'k-', lw=0.25,alpha=0.4)
            

            t_array, forecast, _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast=True, multiplier=m, P_mar=k)
            ax.plot(t_array, forecast, label=f"{m * 100}% of current cow population", lw = 2, alpha = 1)
        ax.plot([],[],'k-', lw=0.5,alpha=0.4, label='model ensemble')
        
        ax.legend()
        #plt.savefig(f"forecast_{j}.jpg")
        plt.show()

def plot_samples2D(a, b, P, samples,b_1, b_2, b_3, tau, p_0, m_0, alpha):
    # plotting



    fig = plt.figure(figsize=[10., 7.])				# open figure
    ax1 = fig.add_subplot(111, projection='3d')		# create 3D axes
    A, B = np.meshgrid(a, b, indexing='ij')
    ax1.plot_surface(A, B, P, rstride=1, cstride=1,cmap=cm.coolwarm, lw = 0.5)	# show surface

    t_calibrate, nitrate_calibrate = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
    v = 1.

    t_numeric, n_numeric, _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha,forecast = True)
    nitrate_cal_int = np.interp(t_numeric, t_calibrate, nitrate_calibrate)

    s = np.array([np.sum((solve_ode(b_1,b_2,b_3, tau, p_0, m_0, alpha, forecast = True)-nitrate_cal_int)**2)/v for a,b in samples])
    p = np.exp(-s/2.)
    p = p/np.max(p)*np.max(P)*1.2
        
    ax1.plot(*samples.T,p,'k.')

    # plotting upkeep
    ax1.set_xlabel('a', fontsize = text_size)
    ax1.set_ylabel('b', fontsize = text_size)
    ax1.set_zlabel('P', fontsize = text_size)
    ax1.set_zlim(0., )
    ax1.view_init(40, 100.)

    # save and show
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






