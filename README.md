# BSA
Bulirsh-Stoer method in plane restricted three body problem

Input files

param.in - parameters of numerical simulation 
beginTime(double) endTime(double) timeStep(double) timeSnapshot(double) kappa(double) 

beginTime - start time; endTime - end point in time; timeStep - time step; timeSnapshot - output step; kappa - planetary accretion radius in units of Hill radius;

stars.in - stars parameters
mass(double) x(double) y(double) vx(double) vy(double) 

particles.in - stars parameters
x(double) y(double) vx(double) vy(double) id(int)

Output files

out_###i_k - particle positions at the time of data output, i - time moment, k - thread number
x(double) y(double) vx(double) vy(double) id(int) dC/C(double) 

remove/remove_k
time(double) x(double) y(double) vx(double) vy(double) id(int) dC/C(double) type(int)

time - time of particle leaving the system
dC/C - relative change in the value of the Jacobi integral
type - 0-accretion, 1-scattering
