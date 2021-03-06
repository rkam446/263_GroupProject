# imports
from math import exp
import numpy as np
from matplotlib import pyplot as plt
from numpy.core.function_base import linspace


def get_n_stock(t, tau, forecast=False, multiplier=1):
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
    #Pull cow and year data from file
    year, cattle = np.genfromtxt("data/nl_cows.txt", delimiter=",", skip_header=1, unpack=True)

    #If function is used to forecast then solve out to 2040
    if forecast:
        final_cattle = cattle[-1] * multiplier
        year = np.append(year, [2040])
        cattle = np.append(cattle, final_cattle)

    #Interpolate cattle values between year values 
    interp_cattle = np.interp(t - tau, year, cattle)
    n_stock = dict(zip(t - tau, interp_cattle))

    return n_stock


def ode_model(t, C, P, n_stock, m_0, t_c, t_mar, P_a, P_mar, b_1, b_2, b_3, tau, alpha, P_s):
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
    #Use either pre 2010 ODE or post 2010 ODE
    b = b_1 * alpha if t - tau > t_c else b_1
        
    P_1 = P_a if t < t_mar else P_a + P_mar
        
    dCdt = (-n_stock * b * (P - P_s) + b_2 * C * (P - P_1 / 2)) / m_0
    dPdt = -b_3 * 2 * P if t < t_mar else -b_3 * (2 * P - P_mar / 2)

    return dCdt, dPdt


def solve_ode(b_1, b_2, b_3, tau, p_0, m_0, alpha, P_s = 0.05, dt=0.2, forecast=False, P_mar=0, multiplier=1, Benchmark = False):
    
    #When used for benchmarking initalise with values used for analytical solution
    if Benchmark:
        steps = int(np.ceil((100/dt)))
        t_array = np.arange(steps + 1) * dt
        n_stock = dict(zip(t_array - tau, [1] * steps))
        P_a = 2
        P_s = 1
    #When used for forecasting solve until 2040.
    elif forecast:
          steps = int(np.ceil((2040 - 1980)/ dt))
          t_array = np.arange(steps + 1) * dt + 1980
          n_stock = get_n_stock(t_array, tau, forecast=True, multiplier=multiplier)
          P_a = 0.1
    #When used for fitting to observed data only solve until 2020
    else: 
         t_nitrate, _ = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
         steps = int(np.ceil((t_nitrate[-1] - t_nitrate[0])/ dt))
         t_array = np.arange(steps + 1) * dt + t_nitrate[0]
         n_stock = get_n_stock(t_array, tau)
         P_a = 0.1

    #Initialise arrat to correct size
    C = 0.*t_array
    P = 0.*t_array

    #Initial conditions of ODE 
    C[0] = 2 if Benchmark else 0.2
    P[0] = p_0


    dCdt_1 = 0.
    dPdt_1 = 0.
    dCdt_2 = 0.
    dPdt_2 = 0.

    for i in range(steps):
        #Calculate gradients at original point
        dCdt_1, dPdt_1 = ode_model(
            t_array[i], C[i], P[i], n_stock[t_array[i] - tau], m_0=m_0, t_c=2010, t_mar=2020, P_a=P_a, P_mar=P_mar, 
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,
            alpha=alpha,
            P_s=P_s
        )
        #Use gradient to find predicted point
        C_1 = C[i] + dt * dCdt_1
        P_1 = P[i] + dt * dPdt_1
        
        #Find gradient of predicted point
        dCdt_2, dPdt_2 = ode_model(
            t_array[i + 1], C_1, P_1, n_stock[t_array[i] - tau], m_0, t_c=2010, t_mar=2020, P_a=P_a, P_mar=P_mar,
            b_1=b_1, 
            b_2=b_2,
            b_3=b_3,
            tau=tau,                                   
            alpha=alpha,
            P_s=P_s
        )
        #Average the gradients of predicted and original point to find corrected gradient
        dCdt = (dCdt_1 + dCdt_2) / 2
        dPdt = (dPdt_1 + dPdt_2) / 2
        
        #Find the next value using the corrected gradient 
        C[i + 1] = C[i] + dt * dCdt
        P[i + 1] = P[i] + dt * dPdt

    return t_array, C, P


def get_nitrate_concentration(t, b_1=1, b_2=1, b_3=1, tau=5, p_0=1, m_0=1e9, alpha=0.5):
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

    #Pull data from file
    t_nitrate, _ = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
    #Solve using data from file 
    t_array, C, _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha)
    #Return interpolated values of the solution
    return np.interp(t_nitrate, t_array, C)


def analytic(t):

    #Initalise solution matrix 
    C = [0.]*t

    #Loop through analytical solution to create array
    for i in range(len(t)):
        C[i] = 1 + np.exp(-t[i]-0.5 * np.exp(-2 * t[i]) + 0.5)

    return C
