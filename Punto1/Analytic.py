import numpy as np
from scipy.optimize import fsolve

# Analytic solution to Sod's Shock Tube problem
# reference: "http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html"

def sod_fun(P):
	# defines function to be used with fsolve
	# initial conditions
	rho_l = 1;
	P_l = 1;
	u_l = 0;

	rho_r = 0.125;
	P_r = 0.1;
	u_r = 0;

	gamma = 1.4;

	mu = np.sqrt( (gamma-1)/(gamma+1) );

	return (P - P_r)*(( ((1 - mu**2)**2)*((rho_r*(P + mu*mu*P_r))**-1) )**(0.5)) - 2*(np.sqrt(gamma)/(gamma - 1))*(1 - P**((gamma - 1)/(2*gamma)));

t=1

# initial conditions
x0 = 2;
rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;

gamma = 1.4;
mu = np.sqrt( (gamma-1)/(gamma+1) );

# speed of sound
c_l = pow( (gamma*P_l/rho_l),0.5);
c_r = pow( (gamma*P_r/rho_r),0.5);

P_post = fsolve(sod_fun,(P_l+P_r)/2.0);
v_post = 2*(np.sqrt(gamma)/(gamma - 1))*(1 - pow(P_post, (gamma - 1)/(2*gamma)));
rho_post = rho_r*(( (P_post/P_r) + mu**2 )/(1 + mu*mu*(P_post/P_r)));
v_shock = v_post*((rho_post/rho_r)/( (rho_post/rho_r) - 1));
rho_middle = (rho_l)*pow((P_post/P_l),1/gamma);

# key Positions
x1 = x0 - c_l*t;
x3 = x0 + v_post*t;
x4 = x0 + v_shock*t;
# determining x2
c_2 = c_l - ((gamma - 1)/2)*v_post;
x2 = x0 + (v_post - c_2)*t;

n_points = 5000;
x_min = 0.0;
x_max = 4.0;

x = np.linspace(x_min,x_max,num=n_points);
X = np.zeros((n_points,1));
rho = np.zeros((n_points,1));
P = np.zeros((n_points,1));
u = np.zeros((n_points,1));
e = np.zeros((n_points,1));

f = open("Analytic.dat", "w")

for index in range(n_points):
    if x[index] < x1:
        #Solution before expansion
        rho[index] = rho_l;
        P[index] = P_l;
        u[index] = u_l;
    elif (x1 <= x[index] and x[index] <= x2):
        #Solution in the expansion fan
        c = mu*mu*((x0 - x[index])/t) + (1 - mu*mu)*c_l; 
        rho[index] = rho_l*pow((c/c_l),2/(gamma - 1));
        P[index] = P_l*pow((rho[index]/rho_l),gamma);
        u[index] = (1 - mu*mu)*( (-(x0-x[index])/t) + c_l);
    elif (x2 <= x[index] and x[index] <= x3):
        #Solution between the expansion and the contact
        rho[index] = rho_middle;
        P[index] = P_post;
        u[index] = v_post;
    elif (x3 <= x[index] and x[index] <= x4):
        #Solution between the contact and the shock
        rho[index] = rho_post;
        P[index] = P_post;
        u[index] = v_post;
    elif x4 < x[index]:
        #Solution after shock
        rho[index] = rho_r;
        P[index] = P_r;
        u[index] = u_r;
    e[index] = P[index]/((gamma - 1)*rho[index]) + u[index]*u[index]/2.0;
    X[index] = index*(x_max-x_min)/n_points+x_min;
    f.write(str(X[index])[1:-1] + "\t" + str(rho[index])[1:-1] + "\t" + str(u[index])[1:-1] + "\t" +
            str(e[index])[1:-1] + "\t" + str(P[index])[1:-1] + "\n")
f.close()
