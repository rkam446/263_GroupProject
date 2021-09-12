import numpy as np
from matplotlib import pyplot as plt
from solve_ode import solve_ode, get_nitrate_concentration, ode_model, get_n_stock, analytic

#Set tolerance for the answer to be within
tol = 2.e-5
t = 0
#Solve using ODE solver
x = get_nitrate_concentration(t, 1,1,1,1,1,1,1)
#Previously calculated answer
ans = 352577.13179164

#Check answer is within tolerance 
check = x[-1] - ans
assert abs(check) < abs(tol)

#Solve using ODE solver
x = get_nitrate_concentration(t, 2,3,5,7,11,13,17)
#Previously calculated answer
ans = 6044894.279571

#Check answer is within tolerance 
check = x[-1] - ans
assert abs(check) < abs(tol)

#Solve using ODE solver
x = solve_ode(1,1,1,1,1,1,1)
#Previously calculated answer
ans = 354730.7819570

#Check answer is within tolerance 
check = x[1][-1] - ans
assert abs(check) < abs(tol)

#Solve using ODE solver
x = solve_ode(2,3,5,7,11,13,17)
#Previously calculated answer
ans = 6187477.006990

#Check answer is within tolerance 
check = x[1][-1] - ans
assert abs(check) < abs(tol)

def Plot_benchmark():
    '''
    No inputs

    No returns

    Prints a graph of the ODE solver vs the analytical solution.
    '''
    t_array, C, *_ = solve_ode(1,1,1,0,1,1,1,0.1,Benchmark=True)
    plt.plot(t_array, analytic(t_array), 'rx', label = 'Analytical Solution')
    plt.plot(t_array, C, 'b', label = 'Numerical Solution')
    plt.xlabel('time')
    plt.legend()
    plt.show()

    t_array, Cstep, *_ = solve_ode(1,1,1,0,1,1,1,dt = 2,Benchmark=True)
    plt.plot(t_array, analytic(t_array), 'rx', label = 'Analytical Solution')
    plt.plot(t_array, Cstep, 'b', label = 'Numerical Solution')
    plt.xlabel('time')
    plt.legend()
    plt.show()