import numpy as np
from scipy.optimize import curve_fit 
import scipy as sc
from matplotlib import pyplot as plt
from solve_ode import get_nitrate_concentration, solve_ode

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

    np.random.seed(10)
    n_samples = 10
    
    values = np.random.normal(loc=paras, scale=perr, size=(n_samples, len(paras)))
    
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
    #plt.savefig("fitted_model.jpg")
    plt.show()
    
    # 4. Plot graphs
    cattle_multipliers = [0, 0.5, 1, 2]

    mar_values = [0, 0.1, 1]
    forecast = list()

    for j, k in enumerate(mar_values): 
        fig, ax = plt.subplots()
        ax.set_xlabel('time (yrs)')
        ax.set_ylabel('nitrate concentration (mg/L)')
        ax.set_title(f"MAR = {k} MPa")

        for i, m in enumerate(cattle_multipliers):
            t_array, forecast, _ = solve_ode(b_1=b_1, b_2=b_2, b_3=b_3, tau=tau, p_0=p_0, m_0=m_0, alpha=alpha, forecast=True, multiplier=m, P_mar=k)
            ax.plot(t_array, forecast, label=f"{m * 100}%")

            for l in range(n_samples):
                t_array, forecast, _ = solve_ode(*values[l, :], forecast=True, multiplier=m, P_mar=k)
                ax.plot(t_array, forecast, alpha=0.2)
        
        ax.legend()
        #plt.savefig(f"forecast_{j}.jpg")
        plt.show()
