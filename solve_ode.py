# imports
from ssl import create_default_context
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
        n_stock : float array
            Array of number of cattle at each time interval.
    """
    year, cattle = np.genfromtxt("data/nl_cows.txt", delimiter=",", skip_header=1, unpack=True)
    interp_cattle = np.interp(t - tau, year, cattle)
    n_stock = dict(zip(t - tau, interp_cattle))
    return n_stock




def ode_model(t, C, P, n_stock, m_0, t_c, t_mar, P_a, P_mar, b_1, b_2, b_3, tau, alpha):
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
        P_a : float
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
    b = b_1 * alpha if t - tau > t_c else b_1
        
    P_1 = P_a if t < t_mar else P_a + P_mar
        
    dCdt = (-n_stock * b * (P - P_a) + b_2 * C * (P - P_1 / 2)) / m_0
    dPdt = -b_3 * 2 * P if t < t_mar else -b_3 * (2 * P - P_mar / 2)

    return dCdt, dPdt
    

def get_nitrate_concentration(t, b_1=1, b_2=1, b_3=1, tau=0.5, p_0=1, m_0=1e9, alpha=0.5):
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

    dt = 0.1

    t_nitrate, _ = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
    t_nitrate = t_nitrate[10:]

    steps = int(np.ceil((t_nitrate[-1] - t_nitrate[0])/ dt))
    t_array = np.arange(steps + 1) * dt + t_nitrate[0]

    p = 0.*t_array
    c = 0.*t_array
    c[0] = 4.2
    p[0] = p_0

    dCdt_1 = 0.
    dPdt_1 = 0.
    dCdt_2 = 0.
    dPdt_2 = 0.

    n_stock = get_n_stock(t_array,tau)

    for i in range(steps):
        dCdt_1, dPdt_1 = ode_model(t_array[i], c[i], p[i], n_stock[t_array[i] - tau], m_0=m_0, t_c=2010, t_mar=2020, P_a=0.05, P_mar=0, 
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,
            alpha=alpha
        )

        c1 = c[i] + dt * dCdt_1
        p1 = p[i] + dt * dPdt_1

        dCdt_2, dPdt_2 = ode_model(t_array[i + 1], c1,p1, n_stock[t_array[i] - tau], m_0, t_c=2010, t_mar=2020, P_a=0.05, P_mar=0,
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,                                   
            alpha=alpha
        )

        dCdt = (dCdt_1 + dCdt_2) / 2
        dPdt = (dPdt_1 + dPdt_2) / 2

        c[i + 1] = c[i] + dt * dCdt
        p[i + 1] = p[i] + dt * dPdt
        
    
    return np.interp(t_nitrate, t_array, c)


def fitted_model(b_1, b_2, b_3, tau, p_0, m_0, alpha):
    dt = 0.1

    t_nitrate, _ = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
    t_nitrate = t_nitrate[10:]

    steps = int(np.ceil((t_nitrate[-1] - t_nitrate[0])/ dt))
    t_array = np.arange(steps + 1) * dt + t_nitrate[0]

    p = 0.*t_array
    c = 0.*t_array
    c[0] = 4.2
    p[0] = p_0

    dCdt_1 = 0.
    dPdt_1 = 0.
    dCdt_2 = 0.
    dPdt_2 = 0.

    n_stock = get_n_stock(t_array,tau)

    for i in range(steps):
        dCdt_1, dPdt_1 = ode_model(t_array[i], c[i], p[i], n_stock[t_array[i] - tau], m_0=m_0, t_c=2010, t_mar=2020, P_a=0.05, P_mar=0, 
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,
            alpha=alpha
        )

        c1 = c[i] + dt * dCdt_1
        p1 = p[i] + dt * dPdt_1

        dCdt_2, dPdt_2 = ode_model(t_array[i + 1], c1,p1, n_stock[t_array[i] - tau], m_0, t_c=2010, t_mar=2020, P_a=0.05, P_mar=0,
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,                                   
            alpha=alpha
        )

        dCdt = (dCdt_1 + dCdt_2) / 2
        dPdt = (dPdt_1 + dPdt_2) / 2

        c[i + 1] = c[i] + dt * dCdt
        p[i + 1] = p[i] + dt * dPdt

    return t_array, c