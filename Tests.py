import numpy as np
from matplotlib import pyplot as plt
from solve_ode import solve_ode, get_nitrate_concentration, ode_model, get_n_stock, analytic

tol = 2.e-5
t = 0
x = get_nitrate_concentration(t, 1,1,1,1,1,1,1)
ans = 352577.13179164

check = x[-1] - ans
assert check < tol

x = get_nitrate_concentration(t, 2,3,5,7,11,13,17)
ans = 6044894.279571

check = x[-1] - ans
assert check < tol


x = solve_ode(1,1,1,1,1,1,1)
ans = 354730.7819570

check = x[1][-1] - ans
assert check < tol

x = solve_ode(2,3,5,7,11,13,17)
ans = 6187477.006990

check = x[1][-1] - ans
assert check < tol

t_array, C, *_ = solve_ode(1,1,1,0,1,1,1,0.1,Benchmark=True)
plt.plot(t_array, analytic(t_array), 'rx')
plt.plot(t_array, C)
plt.show()

t_array, Cstep, *_ = solve_ode(1,1,1,0,1,1,1,dt = 2,Benchmark=True)
plt.plot(t_array, analytic(t_array), 'rx')
plt.plot(t_array, Cstep)
plt.show()