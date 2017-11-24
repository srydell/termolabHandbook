import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import rc

plt.style.use("ggplot")
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

# Real value for heat of vaporisation for nitrogen [kJ/kg]
H = 198.57

# Power in Watt
p1 = [.225, .582, .862, 1.233, 1.692, 5.148, 9.931, 21.235]
p2 = [0.2, .5, .75, 1.25, 1.75, 5, 10, 20]
p3 = [.206, .499, .744, 1.24, 1.744, 4.936, 10.31, 20.235]
p4 = [0.2, .5, .75, 1.25, 1.75, 5, 10, 20]
p5 = [0.2, .5, .75, 1.25, 1.75, 5, 10, 20]
p6 = [0.175, 0.585, 1.212, 1.992, 4.989, 10.54,22.08]
p7 = [0.6314, 1.908, 6.7545, 9.9275, 18.2114, 27.009, 0.33, 4.7158, 0.9316]
p8 = [0.17395, 0.5668, 0.7182, 1.382, 1.84, 4.864, 10.45, 19.798]
p9 = [.1665, .492, 0.7365, 1.2309, 1.6698, 1.7917, 5.15, 13.0758]
p10 = [.11552, .396, .726, 1.312, 1.7934, 5.0388, 10.142, 19.9815]
p11 = [.1771, .396, .75, 1.533, 2.4948, 4.9476, 15.1889, 20.2818]
p12 = [.218, .426, .7695, 1.552, 2.546, 4.967, 15.283, 19.9]
p13 = [.235, .396, .737, 1.491, 2.626, 5.252, 14.848, 20.202]

# Mass difference (mdot) in kg/s * 10^-5
m1 = [2.067, 1.733, 1.8, 1.867, 2.067, 3.2, 5, 9.467]
m2 = [1.84, 2.15, 1.92, 1.97, 2.3, 3.9, 5.7, 10.09]
m3 = [1.73, 1.83, 2, 2.1, 2.33, 3.7, 6.2, 11.2]
m4 = [.96, 1.2, 1.08, 1.4, 1.56, 2.24, 3.64, 6.38]
m5 = [1.66, 1.66, 1.8, 2, 2.86, 4.066, 7.2, 10.22]
m6 = [1.4333, 1.6667,1.7333,1.8667,2.1667,3.5667,6.3667,12]
m7 = [1.73, 2, 4.4, 5.93, 9.93, 14.5, 1.27, 3.33, 1.47]
m8 = [1.3, 1.667, 1.667, 2, 2.23, 3.533, 6.267, 10.8]
m9 = [1.4, 1.5, 1.4, 1.5, 1.2, 1.6, 3.6, 4.3, 7.6]
m10 = [1.53, 1.63, 1.67, 1.8, 2.27, 3.4, 5.93, 10.6]
m11 = [1.3, 2.5, 1.8, 2, 3.5, 3.5, 8.6, 11.1]
m12 = [1.6, 1.6, 1.8, 2.01, 2.47, 3.55, 8.6, 10.9]
m13 = [.917, 1.5, 1.66, 2.166, 2.5, 3.58, 8.75, 11.1]


# Combine mdata into one list
ms = [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13]
mdot = []
for m in ms:
	mdot.extend(m)

# Multiply every value with 10^(-5) to match data taken
mdot = [x * math.pow(10, -5) for x in mdot]

# Combine pdata into one list
ps = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13]
pdot = []
for p in ps:
	pdot.extend(p)

# error ~ 1/P
# error taken as s/J
error = [1/p for p in pdot]

nValues = len(pdot)

# Create 1/P for x-axis and P/mdot for y-axis
# P/mdot ~ enthalpy
# enthalpy taken as kJ/kg
enthalpy = [0]*len(pdot)
for n in range(nValues):
	# Multiply every value with 10^(-3) to get kJ/s
	enthalpy[n] = pdot[n]/mdot[n] * math.pow(10, -3)

# Calculate a best fit
# Add the real value of the enthalpy
# corresponding to 0 error
# enthalpy.extend([H])
# error.extend([0])

# Function to fit against
def func(x, a, b, c):
	return a * np.exp(-b * x) + c

# Fit the data
xdata = error.copy()
xdata.sort()
xdata = np.array(xdata)
ydata = func(xdata, H, 2, 20)

popt, pcov = curve_fit(func, xdata, ydata)

plt.scatter(error, enthalpy, c="#8eb5c6", edgecolors='b', label="Measured data", alpha=0.4, s=40)
plt.plot(xdata, func(xdata, *popt), c="#094f6e", linewidth=2, label=r"$H \cdot e^{-\frac{2}{P}} + 20$")

plt.yticks([0, 100, H, 300], [0, 100, r"$H =\\ {}$".format(H), 300])

plt.xlabel(r"$\frac{1}{P} [\frac{s}{J}$]", fontsize=16)
plt.ylabel(r"$\frac{P}{\dot{m}} [\frac{kJ}{kg}]$", fontsize=16)

plt.xlim(-0.1, xdata[-1] + .1)
plt.ylim(-10, max(enthalpy) + 10)
plt.legend(loc=1, prop={"size": 18})

# Save the plot
plt.savefig('nitrogenEnthalpy.png', bbox_inches='tight')
plt.show()
