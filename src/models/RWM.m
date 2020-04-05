
%-----PHYSICAL PARAMETERS------%
l=334e3;
%specific latent heat [J/kg]
rho=1e-6;
%density [kg/mm^3]
K=0.6e-3;
%thermal conductivity [W/(mm*K)]
c=4.2e3;
%specific heat [J/(kg*K)]
alpha=K/(rho*c); %thermal diffusivity [mm^/s]
beta=c/l; %[1/K]
L=2; %[mm] -Length of domain
t_max=1; %[s] -maximum time
T_0=1; %[degree C] -Temperature at fixed boundary
%-----NUMERICAL VALUES------%
n=1e3; %number of iterations
dx=0.01; %Steplength in x
dt=dx^2/(2*alpha); %Steplength in t
ds=dx/(n*beta); %Increment size for moving boundary when absorbing one walker.
N_x=ceil(L/dx); %Number of points in spatial domain.
N_t=ceil(t_max/dt); %Number of points in time domain.
T=zeros(N_x,N_t); %Matrix representing T(x,t)
s_vector=zeros(1,N_t); %Vector representing moving boundary s(t)
t_j=1;
s_i=1;
s=dx;
while t_j < N_t && s_i < N_x %Loop for all time steps as long as s(t) < L
    %-----Boundary condition for fixed boundary-----%
    T(1,t_j)= T_0*n; %Constant BC
    %T(1,t_j)=(exp(t_j*dt)-1)*n; %Time-dependent exponential BC
    %T(1,t_j)=(sin(t_j*dt))*n; %Time-dependent oscillating BC
    s_vector(t_j)=s;
    for x_i = 1:N_x
        if T(x_i,t_j) < 0 %Handling of temperatures below zero.
            sign=-1;
        else
            sign=1;
        end
        for k=1:sign*T(x_i,t_j) %Move all walkers at position (x_i,t_j)
            p=(-1)^round(rand); %P(+1)=P(-1)=1/2

            if x_i+p > 1 && x_i+p <= s_i && x_i <= N_x %A walker moves if it has not reached the boundaries.
                T(x_i+p,t_j+1) = T(x_i+p,t_j+1) + 1*sign;

            elseif x_i+p == s_i+1 %If a walker reaches the moving boundary the boundary moves ds.
                s=s+ds*sign;
                s_i=floor(s/dx); %Index for s(t) updates when s modulo dx = 0.
            end
        end
    end
    t_j=t_j+1;
end
T=T/n; %Dividing by number of iterations.

%------PLOTS------%
figure(1); hold on %Plots temperature distribution and moving boundary s(t).
[x_matrix,t_matrix] = meshgrid(0:dx:(N_x-1)*dx, 0:dt:(N_t-1)*dt);
mesh(x_matrix',t_matrix',T)
t_vector_end=0:dt:(t_j-1)*dt;
s_t_end=find(s_vector,1,'last')+1;
s_t_vector=s_vector(1:s_t_end);
plot3(s_t_vector,t_vector_end, 0*t_vector_end,'ro')
xlabel('x')
ylabel('t')
zlabel('T')
view([20 40])
set(gca,'FontSize',20)

