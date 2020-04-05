"""
==================
kiln_heatup.py
==================
===============================================================================
A simple script for numerically integrating the time-dependent 1D heat equation
with convective+radiative boundary conditions. This is a very basic simulation
of a kiln during the warmup process.
This is part of the RepKiln project. More information about this code can be
found here:
https://hackaday.io/project/21642-repkiln/log/65537-numerically-solving-the-1d-transient-heat-equation

Version 1.0 Initial Release - August 2017
Copyright 2017 Matt Moses
mmoses152 at gmail dot com
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
===============================================================================
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
from time import sleep

A = 1 # cross sectional area of wall element in m^2
rho = 916.0  # density of wall material in kg / m^3
k = 2.1  # thermal conductivity of wall material in W / (m*K)
c = 2.097 * 1000  # specific heat capacity in J / (kg*K)
sigma = 5.6704E-08  # Stefan-Boltzmann constant in W * m^-2 * K^-4
h = 0.0  # convective heat transfer coefficient in W / (m^2 * K)
Q_dot_in = -300.0 # heater power in watts

T_initial = 263.0 # initial temperature in Kelvin
T_inf = 273.0 # ambient temperature in Kelvin

L = 0.01 # thickness of the entire wall in meters
N = 2 # number of discrete wall segments
dx = L/N # length of each wall segment in meters

total_time = 5*60.0 # total duration of simulation in seconds
nsteps = 500 # number of timesteps
dt = total_time/nsteps # duration of timestep in seconds

# The size of this nondimensional factor gives a rough idea
# of the stability of the crude numerical integration we are using.
# If this is too big there will be problems.
simfac = (k*dt) / (c*rho*dx*dx)
print(simfac)

# this is the factor by which to multiply Q_dot_in and Q_dot_out
heatfac = dx / (k*A)

# initialize volume element coordinates and time samples
x = np.linspace(0, dx*(N-1), N)
timesamps = np.linspace(0, dt*nsteps, nsteps+1)

# create a 2D mesh grid for plotting
# a note on using meshgrid from
# https://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html
#
# "In the 2-D case with inputs of length M and N, the outputs are of
# shape (N, M) for 'xy' indexing and (M, N) for 'ij' indexing."
X, TIME = np.meshgrid(x, timesamps, indexing='ij')

# initialize a big 2D array to store temperature values
T = np.zeros((X.shape))

# set the initial temperature profile of the wall
for i in range(len(x)):
   T[i, 0] = T_initial + i

# Debug stuff, uncomment if needed:
# print X.shape
# print timesamps
# print x
# print TIME
# print X
# print T


#T_out = T[len(x)-1, 0]
#Q_dot_out = sigma * A * (pow(T_out,4) - pow(T_inf,4)) + h * A * (T_out - T_inf)

#print Q_dot_in
#print Q_dot_out

for j in tqdm(range(len(timesamps)-1)):

   # get the outside wall temperature and heat flow at current time
   T_out = T[len(x)-1, j]
   Q_dot_out = sigma * A * (pow(T_out,4) - pow(T_inf,4)) + h * A * (T_out - T_inf)
   # now compute temperature at the outside boundary for the next time step
   T[len(x)-1, j+1] = T_out + simfac * (T[len(x)-2, j] - T_out - heatfac * Q_dot_out)
   # and now compute temperature at the inside boundary for the next time step
   T[0, j+1] = T[0,j] + simfac * (T[1,j] - T[0,j] + heatfac * Q_dot_in)
   # now loop through the interior elements to get their temp for the next time
   for i in range(len(x)-2):
      T[i+1,j+1] = T[i+1,j] + simfac * (T[i,j] - 2*T[i+1,j] + T[i+2,j])

print(T[:,- 1].mean())
# this plots the temperature vs time data as a surface
fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
foo = ax.plot_surface(X*100, TIME/3600, T)


ax.set_xlabel('x (cm)')
ax.set_ylabel('time (hr)')
ax.set_zlabel('temperature (K)')
plt.show()