/*
 * A ___shadow rate___ New Keynesian Model for __Monetary Economies__
 * Considering a ____inflation shock______ 
 */


var  y,c,ce,m,k,i, s, pi,x,b,a, tec_shock, pref_shock, inf_shock, beta_t; 
varexo eps_a, eps_beta, eps_pi;

parameters beta, gamma, sigma, omega, psi, chi, delta, alpha, eta_c, 
      eta_r, phi_s, phi_y, phi_pi, X, LV,rp,rho_a, rho_beta, s_ss,rr_ss;

beta  = 0.99; //Discount factor
gamma = 0.98; //Discount factor of entrepreneurs
sigma = 1 ; // Intertemporal elasticity of substitution
omega = 0.47; //Output elasticity of (real)marginal cost
psi   = 0.5; //Degree of price stickiness

/*
 * This parameter changes depending if monetary or cashless
 */
chi   = 0.02; // Degree of non-separability of utility function
delta = 0.03; // Capital depreciation rate
alpha = 0.3 ; // Capital share in production


eta_c = 1   ; // Consumption elasticity of money demand
eta_r = 28 ; // Interest rate semi-elasticity of money demand


phi_s = 0.73; // Interest rate persistance
phi_y = 0.13; // Interest rate response to output
phi_pi= 1.27; // Interest rate response to inflation

X     = 1.05; // steady-state gross markup
LV    = 0.89; // loan-to-value ratio
rp = 1.009; //Constant risk premium


rho_a = 0.90; // Autocorrelation of technology shock
rho_beta = 0.80 ; //Autocorrelation of preference shock

s_ss = 1/(beta*rp);
rr_ss = 1/beta;

model;
%Most of the steady-state parameters depend of each other, and as all are relative to Output, that variable is omitted.
#lambda = (1-psi)*(1-beta*psi)/psi;
#I = (1/X) *(gamma*alpha*delta)/(1-beta*LV+gamma*(LV-1+delta));
#B = LV*I*beta/delta;
#C = ((X-alpha)/X)+B*(1-beta)/beta;
#CE = 1-C-I; 


[name = 'Aggregate Demand']
y = C*c+CE*ce+I*i;
[name = 'Taylor rule']
s = phi_s*s(-1)+(1-phi_s)*(phi_pi*pi+phi_y*y);

[name = 'IS Curve']
c = c(+1)-sigma*(s+beta_t(+1)-pi(+1)+chi*(m(+1)-m));
[name = 'LM Curve']
m = eta_c*c-eta_r*s;

[name = 'Entrepreneur restriction']
CE*ce = (alpha*(y-x)/X)+B*b-I*i-(s(-1)+b(-1)-pi)*B/beta;
[name = 'Lending constraint']
b = pi(+1)+k-s;
[name = 'Entrepreneur Euler']
0 = (1-LV*beta)*(ce-ce(+1))+(LV*beta)*(pi(+1)-s)+(alpha*gamma*delta/(X*I))*(y(+1)-x(+1)-k);

[name = 'Production function']
y = (1+omega)*(a + alpha*k(-1))/(alpha+omega)-(1-alpha)*(x+c/sigma-chi*m)/(alpha+omega);
[name = 'Phillips Curve']
pi = beta*pi(+1)-lambda*x+eps_pi;
[name = 'Law of motion']
k = (1-delta)*k(-1)+delta*i;

[name = 'Technology shock']
a = rho_a*a(-1)+eps_a;
[name = 'Preference shock']
beta_t = rho_beta*beta_t(-1)+eps_beta;

tec_shock = eps_a;
pref_shock = eps_beta;
inf_shock = eps_pi;

end;


initval;
s = 0;
pi = 0;
y = 0;
c = 0;
ce = 0;
m = 0;
k = 0;
i = 0;
x = 0;
b = 0;
end;


check;
steady;

% Define shocks
shocks;
var eps_beta;
periods 1:15;
values (0.05/15);

var eps_pi;
periods 10;
values -0.002;
end;

perfect_foresight_setup(periods=40);
perfect_foresight_solver(lmmcp);