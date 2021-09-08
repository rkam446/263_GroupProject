# imports
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit 

def get_n_stock(t, tau):
    """
        Get an array of cattle numbers for every point on the time interval t.

        Parameters
        ----------
        t : float array
            Time intervals.

        Returns
        -------
        n_wtock : float array
            Array of number of cattle at each time interval.
    """
    n = np.genfromtxt("data/nl_cows.txt", delimiter=",", skip_header=1)
    cows = np.interp(t - tau, n[:, 0], n[:, 1])
    n_stock = dict(zip(t - tau, cows))
    return n_stock




def ode_model(t, C, P, n_stock, m_0, t_c, t_mar, P_surface, P_a, P_mar, b_1, b_2, b_3, tau, alpha):
    """
        Evaluate the rate of change of nitrate concentration and pressure inside the CV at time t.
        
        Parameters
        ----------
        t : float 
            Time of evaluation.
        C : float
            Nitrate concentration in the CV at t.
        P : float
            Pressure inside the CV at t.
        m_0 : float
            Initial mass inside CV.
        t_c : float
            Time at which carbon sink was installed.
        t_mar : float
            Time at which MAR starts.
        P_surface : float
            Pressure at the surface boundary.
        P_a : float
            Base pressure at the CV boundary.
        P_mar : float
            Pressure increase at the CV boundary due to MAR.
        b_1, b_2, b_3, tau, alpha : float
            Lump parameters/constants.
        
        Returns
        -------
        dCdt : float
            Rate of change of C at t.
        dPdt : float
            Rate of change of P at t.
    """
    b = b_1 if t - tau > t_c else alpha * b_1
        
    P_1 = P_a if t < t_mar else P_a + P_mar
        
    dCdt = (-n_stock[t - tau] * b * (P - P_surface) + b_2 * C * (P - P_1 / 2)) / m_0
    dPdt = -b_3 * 2 * P if t < t_mar else -b_3 * (2 * P - P_mar / 2)

    return dCdt, dPdt
    

def get_nitrate_concentration(t, b_1, b_2, b_3, tau, p_0, m_0, alpha):
    ''' Get numeric estimation of the nitrate concentration for a certain year
        Parameters:
        -----------
        t : float
            Time at which to evaluate solution.
        b_1, b_2, b_3, tau, alpha : float
            Parameters/constants.
        Returns:
        --------
        x : float
            Estimated nitrate concentration for the year.
    '''
    t0 = 1980
    t1 = t
    dt = 1
    x0 = [0.2, p_0] # [initial nitrate concentration, initial pressure]

    steps = int(np.ceil((t1 - t0)/ dt))
    t_array = np.arange(steps + 1) * dt
    x = 0.*t_array
    x[0] = x0

    n = get_n_stock(t_array)

    for i in range(steps):
        f0 = ode_model(t_array[i], *x[i], n, m_0, t_c=2010, t_mar=2020, P_surface=0.1, P_mar=0, 
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,
            alpha=alpha
        )

        x1 = x[i] + dt * f0
        f1 = ode_model(t_array[i + 1], x1, n, m_0, t_c=2010, t_mar=2020, P_surface=0.1, P_mar=0, 
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,                                   
            alpha=alpha
        )

        f2 = (f0 + f1) / 2
        x[i + 1] = x[i] + dt * f2
        
    
    return t_array, x[-1]


