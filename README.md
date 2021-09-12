# 263 Team Nine
***This is the Team 9 GitHub repositry for the Edendale Aquifer model.***
The following model is used to model the level of nitrate concentration within the Southland Aquifer.
It uses an Improved Euler solver and the curve_fit function from the scipy library to calibrate unkown model parameters to fit the dataset provided.
It is used to make forecast out to the year 2040 using a variety of what-if scenarios to model the effectiveness of various nitrate leaching mitigating initiatives.
When the *main.py* script is run it calls functions from *solve_ode.py, uncertainty.py, and Tests.py* to create and generate plots of interest.
The *data* folder contains the dataset files used to calibrate the model and are used by some of the functions.
