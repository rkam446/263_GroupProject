import numpy as np
from numpy.core.function_base import linspace
import scipy
import solve_ode
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as curve_fit

if __name__ == "__main__":
    # 1. Get the nitrate concentration data 
    # 2. Get optimal parameters using curve_fit()
    # 3. Solve ODE numerically using optimal parameters
    # 4. Plot numeric solution and actual data
    t = linspace(start=1980, endpoint=2020)
    cows = solve_ode.get_n_stock(t, tau)
    nitrate_conc = solve_ode.get_nitrate_concentration(t= ,b_1= , b_2= , b_3= , tau= , alpha= )

    fit = curve_fit(solve_ode.ode_model, t, nitrate_conc)



    plt.plot(t, nitrate_conc, 'r-', label = 'Model of Nitrate Concentration in Southland Aquifer')
    plt.plot(t, cows)
    plt.show()

    pass
