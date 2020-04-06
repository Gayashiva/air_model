import pandas as pd
import numpy as np
import os
import math
from tqdm import tqdm
import time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from src.data.config import site, dates, option, folders, fountain, surface

#Stefan_problem_RW.m uses a random walk method to solve the heat equation
#with inhomogenous BC on 0 < x < s(t), t>0, where s(t) is a moving boundary.
#BC1: T(0,t) = f(t), t>0 (arbitrary function)
#BC2: T(s(t),t)=0 (melting temperature)
#IC: T(x,0)=0
#Stefan condition: l*rho*ds/dt=K*dT(s(t),t)/dx
#s(0)=0 (here set to dx approx. 0)

#-----PHYSICAL PARAMETERS------#
l=334e3  #specific latent heat [J/kg]

rho=1e-6  #density [kg/mm^3]
K=0.6e-3  #thermal conductivity [W/(mm*K)]
c=4.2e3  #specific heat [J/(kg*K)]

# rho=0.916e-6  #density [kg/mm^3]
# K=2.22e-3  #thermal conductivity [W/(mm*K)]
# c=2.1e3  #specific heat [J/(kg*K)]

alpha=K/(rho*c)  #thermal diffusivity [mm^/s]
beta=c/l  #[1/K]
L=30  #[mm] -Length of domain
t_max= 5 * 60  #[s] -maximum time

#-----NUMERICAL VALUES------ #
n=1e2   #number of iterations
dx=0.1   #Steplength in x
dt=dx**2/(2*alpha)   #Steplength in t
ds=dx/(n*beta)   #Increment size for moving boundary when absorbing one walker.
N_x= math.ceil(L/dx)   #Number of points in spatial domain.
N_t= math.ceil(t_max/dt)   #Number of points in time domain.
T= np.zeros((N_x,N_t))   #Matrix representing T(x,t)
s_vector= np.zeros((N_t))   #Vector representing moving boundary s(t)
t_j=0
s_i=0
s=dx

# Q = 100e-6 # W/mm2
# T_0=  Q * dx/K #[degree C] -Temperature at fixed boundary
T_0 = 0.1
print(T_0)

pbar = tqdm(total=N_t)
while t_j < N_t-1 and s_i < N_x-1:

    # -----Boundary condition for fixed boundary-----#
    T[0, t_j] = T_0 * n
    # T(1,t_j)=(exp(t_j*dt)-1)*n; #Time-dependent exponential BC
    # T(1,t_j)=(sin(t_j*dt))*n; #Time-dependent oscillating BC
    s_vector[t_j] = s

    for x_i in np.arange(0, N_x):
        if T[x_i, t_j] < 0:
            sign = - 1
        else:
            sign = 1

        for k in np.arange(0, sign*T[x_i, t_j]):
            p = (- 1) ** round(np.random.uniform(0, 1))

            if x_i + p > 0 and x_i + p <= s_i and x_i <= N_x-1:

                T[x_i + p, t_j + 1] = T[x_i + p, t_j + 1] + sign

            else:

                if x_i + p == s_i + 1:
                    s = s + ds*sign

                    s_i = math.floor(s / dx)

    pbar.update(1)
    t_j = t_j + 1

pbar.close()
T = T / n

# initialize volume element coordinates and time samples
x = np.linspace(0, dx * N_x, N_x)
timesamps = np.linspace(0, dt * N_t, N_t)

# "In the 2-D case with inputs of length M and N, the outputs are of
# shape (N, M) for 'xy' indexing and (M, N) for 'ij' indexing."
X, TIME = np.meshgrid(x, timesamps, indexing='ij')

t_vector_end = np.arange(0, (t_j-1)*dt, dt)

s_t_end = np.nonzero(s_vector)[-1][-1]

s_t_vector = s_vector[0: s_t_end]

print(s_t_vector[-1])

# ------PLOTS------#
fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
foo = ax.plot_surface(X, TIME, T)
ax.scatter(s_t_vector, t_vector_end, 0,color='red')
ax.set_xlabel('x (mm)')
ax.set_ylabel('time (s)')
ax.set_zlabel('temperature (C)')
# plt.savefig(
#         os.path.join(folders["output_folder"], "Stefan.jpg"),
#         bbox_inches="tight",
#         dpi=300,
#     )
plt.show()

# plt.figure(1)
# hold('on')
# x_matrix, t_matrix = np.meshgrid(np.arange(0, np.dot((N_x - 1), dx), dx), np.arange(0, np.dot((N_t - 1), dt), dt), nargout=2)
# # RWM.m:58
# mesh(x_matrix.T, t_matrix.T, T)
# t_vector_end = np.arange(0, np.dot((t_j - 1), dt), dt)
# # # RWM.m:60
# s_t_end = find(s_vector, 1, 'last') + 1
# # RWM.m:61
# s_t_vector = s_vector(np.arange(1, s_t_end))
# # RWM.m:62
# plot3(s_t_vector, t_vector_end, np.dot(0, t_vector_end), 'ro')
# xlabel('x')
# ylabel('t')
# zlabel('T')
# view(concat([20, 40]))
# set(gca, 'FontSize', 20)