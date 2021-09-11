import numpy as np
from solve_ode import fitted_model, get_nitrate_concentration, ode_model, get_n_stock

tol = 2.e-5
t = 0
x = get_nitrate_concentration(t, 1,1,1,1,1,1,1)
ans = 336775.01081946

check = x[-1] - ans
assert check < tol

x = get_nitrate_concentration(t, 2,3,5,7,11,13,17)
ans = 227248.43362524

check = x[-1] - ans
assert check < tol

tol = 2.e-5
t = 0
x = fitted_model(1,1,1,1,1,1,1)
ans = 338495.4669559

check = x[1][-1] - ans
assert check < tol

x = fitted_model(2,3,5,7,11,13,17)
ans = 233258.52140975

check = x[1][-1] - ans
assert check < tol
