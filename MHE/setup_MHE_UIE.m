function [solver, param_out] = setup_MHE_UIE(  config, solver_options)

%% MHE use current measurment
import  casadi.*

n_y=config.n_y;
n_x=config.n_x;
n_d = config.n_d;
n_var = n_y+n_x+n_d;
n_par = config.n_y_m+config.nx+config.n_u + config.n_factor+ config.n_Nk;
T_tr = config.T;

% optimization variables
z=MX.sym('xw', n_var);
x = reshape(z(1:n_x), [config.nx, config.N+1]);   % est.states
d = reshape(z(n_x+1: n_d+n_x), [config.nd, config.N+1]); % est. disturbance/noices
y = reshape(z(n_x+n_d+1:end), [config.ny, config.N+1]); % est.output
 
% parameter for measurment 
par = MX.sym('par', n_par);
x_int = par(1:config.nx);
y_m = reshape( par(config.nx +1:config.n_y_m+config.nx), [config.ny_m, config.N+1]);
Nk = par(config.n_y_m+config.nx+1:config.n_y_m+config.nx+config.n_Nk);
factor = par(config.n_y_m+config.nx+config.n_Nk+1: config.n_y_m+config.nx+config.n_Nk+config.n_factor);
%u  = reshape(par(config.n_y_m+config.nx+config.n_Nk+1:end), [config.nu, config.N]);
u   = par(config.n_y_m+config.nx+config.n_Nk+config.n_factor+1:end);



obj = MX(0);
ceq = [];
cineq = [];


%% cost and dynamics
for i = 1:config.N

%% dynamic
 ceq =[ceq; factor(i)*x(:,i+1)-factor(i)*config.f_dyn(x(:,i), d(:,i), y_m(:,i), u)];
 ceq = [ceq; factor(i)*y(:,i)-factor(i)*config.h_dyn(x(:,i),d(:,i),y_m(:,i))];
  
%% state constraints    
cineq = [cineq; factor(i)*T_tr(x(:,i),y_m(:,i),d(:,i))-config.xlim(:,2)];
cineq = [cineq; factor(i)*config.xlim(:,1)-factor(i)*T_tr(x(:,i),y_m(:,i),d(:,i))];
%% disturbance and noises constraints
 cineq = [cineq; d(:,i)-config.dlim(:,2)];
 cineq = [cineq; config.dlim(:,1)-d(:,i)];
%% debug
%cineq = [cineq; d(2:3,i)-config.dlim(2:3,2)];
%cineq = [cineq; config.dlim(2:3,1)-d(2:3,i)];
% ceq = [ceq; d(1,i)];  % debug

% ceq = [ceq; d(:,i)];  % debug
%% output constraints
cineq = [cineq; factor(i)*y(:,i)-config.ylim(:,2)];
cineq = [cineq; factor(i)*config.ylim(:,1)-factor(i)*y(:,i)];


% cost
obj = ((config.alpha-1)*factor(i)+1)*obj +  config.alpha*factor(i)*config.cost_dy(y_m(:,i), y(:,i), d(:, i) );


end
%% terminal disturbance and noises constraints
cineq = [cineq; d(:,Nk+1)-config.dlim(:,2)];  % add future disturbance constraints
cineq = [cineq; config.dlim(:,1)-d(:,Nk+1)];

%% terminal output model  
ceq = [ceq; y(:,Nk+1)-config.h_dyn(x(:,Nk+1),d(:,Nk+1),y_m(:,Nk+1))]; % add future output function


%% add terminal output constraint 
cineq = [cineq; y(:,Nk+1)-config.ylim(:,2)];  % add future output constraints
cineq = [cineq; config.ylim(:,1)-y(:,Nk+1)];


%% add terminal stage cost
obj = obj+ config.cost_dy(y_m(:,1+Nk), y(:,Nk+1), d(:, 1+Nk) );  % add future measurement into tht cost


%% terminal state constraint
cineq = [cineq; x(:,Nk+1)-config.xlim_t(:,2)];
cineq = [cineq; config.xlim_t(:,1)-x(:,Nk+1)];
% add prior cost
obj = obj + config.cost_init(Nk,x_int, x(:,1));



param_out.n_eq = length(ceq);
param_out.n_ineq = length(cineq);
param_out.n_var = n_var;

con=[ceq; cineq];
nlp = struct('x', z, 'p', par, 'f', obj, 'g', con );
solver = nlpsol('solver', 'ipopt', nlp, solver_options);
end



