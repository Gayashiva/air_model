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
alpha=K/(rho*c)  #thermal diffusivity [mm^/s]
beta=c/l  #[1/K]
L=2  #[mm] -Length of domain
t_max=1  #[s] -maximum time
T_0=1  #[degree C] -Temperature at fixed boundary

#-----NUMERICAL VALUES------ #
n=1e3   #number of iterations
dx=0.01   #Steplength in x
dt=dx**2/(2*alpha)   #Steplength in t
ds=dx/(n*beta)   #Increment size for moving boundary when absorbing one walker.
N_x= ceil(L/dx)   #Number of points in spatial domain.
N_t= ceil(t_max/dt)   #Number of points in time domain.
T= np.zeros(N_x,N_t)   #Matrix representing T(x,t)
s_vector= np.zeros(1,N_t)   #Vector representing moving boundary s(t)
t_j=1 
s_i=1 
s=dx 

while (t_j < N_t) and (s_i < N_x) :  # Loop for all time steps as long as s(t) < L

     # -----Boundary condition for fixed boundary - ----  #
    T(1, t_j) = T_0 * n  # Constant BC
    # T(1, t_j) = (exp(t_j * dt) - 1) * n;  # Time - dependent exponentialBC
    # T(1, t_j) = (sin(t_j * dt)) * n;  # Time - dependent oscillating BC

    s_vector(t_j) = s

    for x_i in range(1,N_x+1):  #Handling of temperatures below zero.

        if T(x_i,t_j) < 0:
            sign = -1
        else:
            sign =1

        for k in range(1, sign*T(x_i,t_j))  #Move all walkers at position (x_i,t_j)
            p = (-1) ^ round(rand)  # P(+1) = P(-1) = 1 / 2

            if (x_i + p > 1) and (x_i + p) <= s_i and (x_i <= N_x):  # A walker moves if it has not reached the boundaries.
                T(x_i + p, t_j + 1) = T(x_i + p, t_j + 1) + 1 * sign

            else if:
            x_i + p == s_i + 1  # If a walker reaches the moving boundary the boundary moves ds.
