import numpy as np
from matplotlib import pyplot as plt

textsize = 16

t_nitrate, nitrate = np.genfromtxt("data/nl_n.csv", delimiter=",", skip_header=1, unpack=True)
t_cattle, cattle = np.genfromtxt("data/nl_cows.txt", delimiter=",", skip_header=1, unpack=True)

fig, ax1 = plt.subplots()

ax1.set_xlabel("Year")
ax1.set_ylabel("Nitrate Concentration (mg/L)")
ax1.scatter(t_nitrate, nitrate, label="Nitrate", c="red")
ax1.scatter([], [], label="Cattle", c="black")
v = 1.
ax1.errorbar(t_nitrate,nitrate,yerr=v,fmt='ro', label='data')
ax2 = ax1.twinx()
ax2.set_ylabel("Cattle Numbers")
ax2.scatter(t_cattle, cattle, label="Cattle", c="black")
ax1.legend()
plt.show()