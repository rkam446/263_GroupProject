import numpy as np
import scipy as sc
from matplotlib import pyplot as plt
from solve_ode import fitted_model, get_nitrate_concentration, ode_model, get_n_stock

if __name__ == "__main__":
    # 1. Get the nitrate concentration data
    #t = ? #final time of model
    #C = ?
    #P = ?
    #m_0 = ?
    #t_c = ?
    #P_surface = ?
    #P_a = ?
    #P_mar = ?
    
    t_calibrate, nitrate_calibrate = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
    t_calibrate = t_calibrate[10:]
    nitrate_calibrate = nitrate_calibrate[10:]
    # 2. Get optimal parameters using curve_fit()
    
    #n_stock = get_n_stock(t_calibrate)
    #get stock numbers for t_calibrate

    paras, _ = sc.optimize.curve_fit(get_nitrate_concentration, xdata=t_calibrate, ydata=nitrate_calibrate)
    [b_1, b_2, b_3, tau, p_0, m_0, alpha] = paras
    #? is intial guess
    
    # 3. Solve ODE numerically using optimal parameters
    t_array, nitrate_conc_num = fitted_model(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha)
    
    # 4. Plot numeric solution and actual data

    #plot data predicted by model
    plt.plot(t_array, nitrate_conc_num, label="Modelled")
    #plot measured data
    plt.scatter(t_calibrate, nitrate_calibrate, label="Experimental")

    #label graph
    plt.xlabel('time (yrs)')
    plt.ylabel('nitrate concentration (mg/L)')
    plt.legend()
    plt.title('Experimental vs Modelled Nitrate Concentration for Southland Aquifer')

    plt.show()
    
