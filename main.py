import numpy as np
from numpy.core.function_base import linspace
from scipy.optimize import curve_fit 
import scipy as sc
from matplotlib import pyplot as plt
from solve_ode import get_nitrate_concentration, solve_ode, analytic
from uncertainty import grid_search, construct_samples, model_ensemble
from Tests import Plot_benchmark

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

    ps = np.random.multivariate_normal(paras[0:2], pcov[0:2,0:2], 10)   # samples from posterior
    a,b, posterior = grid_search(b_1, b_2, b_3, tau, p_0, m_0, alpha)
    N = 10
    samples = construct_samples(a, b, posterior, N, b_1, b_2, b_3, tau, p_0, m_0, alpha)

    # Produce Model Predictions with uncertainty.
    model_ensemble(samples,b_1, b_2, b_3, tau, p_0, m_0, alpha)

    # 3. Solve ODE numerically using optimal parameters
    t_array, n_numeric, *_ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha)
    t_array, no_carbon, *_ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=1)
    #Plot solutions on graphs
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time (yrs)')
    ax1.set_ylabel('nitrate concentration (mg/L)')
    ax1.set_title("Fitted Model")
    ax1.plot(t_array, n_numeric, "c-", label = 'With Carbon Sink')
    ax1.plot(t_array, no_carbon, 'b-', label = 'Without Carbon Sink')
    ax1.scatter(t_calibrate, nitrate_calibrate, c="red", label = 'Recorded Nitrate Concentration')
    ax2 = ax1.twinx()
    ax2.set_ylabel("cattle number")
    ax2.scatter(year, cattle, c="black")
    ax1.legend()
    plt.savefig("fitted_model.jpg")
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time (yrs)')
    ax1.set_ylabel('nitrate concentration (mg/L)')
    ax1.set_title("Residuals")
    #Calculate residuals by taking the difference between model and observations
    ax1.scatter(t_calibrate, get_nitrate_concentration(t_array, b_1, b_2, b_3, tau, p_0, m_0, alpha) - nitrate_calibrate)
    ax1.plot(t_calibrate, [0]*len(t_calibrate), 'r-')
    ax1.legend()
    plt.savefig("fitted_model.jpg")
    plt.show()


    Plot_benchmark()
    



    