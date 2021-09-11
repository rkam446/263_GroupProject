import numpy as np
from numpy.core.function_base import linspace
from scipy.optimize import curve_fit 
import scipy as sc
from matplotlib import pyplot as plt
from solve_ode import get_nitrate_concentration, solve_ode, analytic

if __name__ == "__main__":
    # 1. Get measurement data
    t_calibrate, nitrate_calibrate = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
    year, cattle = np.genfromtxt("data/nl_cows.txt", delimiter=",", skip_header=1, unpack=True)
    # t_calibrate = t_calibrate[10:]
    # nitrate_calibrate = nitrate_calibrate[10:]

    # 2. Get optimal parameters using curve_fit()
    #p_0 =[0.001, 24.897, 2.091, 5, -0.650, 14.35, 0.19]
    paras, pcov = sc.optimize.curve_fit(get_nitrate_concentration, xdata=t_calibrate, ydata=nitrate_calibrate)
    [b_1, b_2, b_3, tau, p_0, m_0, alpha] = paras
    perr = np.sqrt(np.diag(pcov))

    fig, ax = plt.subplots()
    ps = np.random.multivariate_normal(paras[0:2], pcov[0:2,0:2], 10)   # samples from posterior
    

    


    # 3. Solve ODE numerically using optimal parameters
    t_array, n_numeric, _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha)
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time (yrs)')
    ax1.set_ylabel('nitrate concentration (mg/L)')
    ax1.set_title("Fitted Model")
    ax1.plot(t_array, n_numeric, "c-")
    ax1.scatter(t_calibrate, nitrate_calibrate, c="red")
    ax2 = ax1.twinx()
    ax2.set_ylabel("cattle number")
    ax2.scatter(year, cattle, c="black")
    plt.savefig("fitted_model.jpg")
    plt.show()

    # 4. Plot graphs
    cattle_multipliers = [0, 0.5, 1, 2]

    mar_values = [0, 0.1, 1]
    forecasts = [[0]*len(cattle_multipliers)]*len(mar_values)
    for j, k in enumerate(mar_values):    
        fig, ax = plt.subplots()
        ax.set_xlabel('time (yrs)')
        ax.set_ylabel('nitrate concentration (mg/L)')
        ax.set_title(f"MAR = {k} MPa")
        for i, m in enumerate(cattle_multipliers):
            t_array, forecasts[j][i], _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast=True, multiplier=m, P_mar=k)
            ax.plot(t_array, forecasts[j][i], label=f"{m * 100}%")
        
        ax.legend()    
        plt.savefig(f"forecast_{j}.jpg")
        plt.show()



    t_array, C, *_ = solve_ode(1,1,1,0,1,1,1,0.1,Benchmark=True)



    plt.plot(t_array, analytic(t_array), 'rx')
    plt.plot(t_array, C)
    plt.show()