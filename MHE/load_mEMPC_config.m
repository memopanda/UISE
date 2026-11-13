
config.Time = 2; %[h] horizon time in h
config.N=config.Time*60*60/plant.dt;  % horizon steps
config.Time_Sim=20; %[h]
config.N_sim = config.Time_Sim*60*60/plant.dt; %[h]-->steps

% config.wx_init = [0.18; 0.0015 ; 30];  % Xd; Xc; X_T (biomass(>0.08, <0.15), CO2(<0.00275), Temperature[5 35])
% config.wx_init = [0.18; 0.0015 ; 15];  % Xd; Xc; X_T (biomass(>0.08, <0.15), CO2(<0.00275), Temperature[5 35])
if plant.nw==1
config.wx_init = [0.1; 0.0005 ; 15];  % Xd; Xc; X_T (biomass(>0.07), CO2(<0.00275), Temperature[5 15])
else
%  config.wx_init = [0.12; 0.12; 0.0025 ; 15];    % old
%   config.wx_init = [0.12; 0.12; 0.0001 ; 8];  % new
  config.wx_init = [0.12; 0.12; 0.0005 ; 6];  % new 2
  config.wx_init = [0.12; 0.12; 0.0005 ; 5.05];  % new 2-->paper
%    config.wx_init = [0.12; 0.12; 0.0005 ; 6.5];  % for new constraints of x_T 

%    config.wx_init = [0.12; 0.12; 0.0001 ; 5.1]; 
end
config.cost_w = 10; %120
config.cost_w_mEMPC = 10; %180
% config.cost = @(w,u, Q) (-Q*(w)'*ones(plant.nw,plant.nw)*(w)+u'*diag([0.01 0.5 2])*u);  % u_v, u_h, u_c  

config.cost = @(w,u, Q)(-Q*(w)'*ones(plant.nw,plant.nw)*(w)+u'*diag([0.01 0.5 0.5])*u);  % u_v, u_h, u_c --> paper

options = struct;
options.ipopt.max_iter = 800;
options.ipopt.print_level = 1;
options.print_time = 0;


config.n_v = config.N*plant.nu;  % no terminal constraint
config.n_v_t= plant.nu+1;  %  terminal constraint nu+1 (n_h1, u_h2, u_v u_c)

config.n_x =(config.N+1)*plant.n;
config.n_e = 0; 
config.xlim  = [0 0.00275; 5 35];  % old: [0 0.00275; 5 35] 
% config.xlim  = [0 0.00275; 13 35];  % old: [0 0.00275; 5 35] 
%debug
% config.xlim  = [0 300; -5 400];
if plant.nu==3
  config.ulim  = [0 1; 0 1; 0 1];  
else
config.ulim  = [0 1; 0 1];
end

if mEMPC_controller
config.n_w = (config.Time*60*60/plant.DT+1)*plant.nw;
else
  config.n_w = (config.N+1)*plant.nw;
end


% for time-vary outside temprature [2~10 Grad]
% T_ver= [ 3 3 3 4 5 5.5 6 6 6.7 7 8 8 8 8  8 8.5 8.5 9 9  9 9.3 9.5 ];
T_ver= [ 3.8 2.5 2 3 4 5.5 6 6 6.7 7 8 8 7.1  6.5  6.5  5.5  5.8  4.8  5 4.6 4.3 3.8 ]; %->paper
% T_ver= [ 3.8 3.3 3 4 5 5.5 6 6 6.7 8 8 9 10  10  9.7  9.5  9  8.7  8  8  7.3 7 ]; %
% T_ver= [ 3.8 3.3 3 4 5 5.5 6 6 6.7 5 5 4.4 4.1  3.8  3.2  3.7  3.1  2.8  2.5 2.3 3.1 3.2 ];
% T_ver = 10*ones(8);
nT=config.N/config.Time;
T_len = (config.Time_Sim+config.Time)*nT;
T_traj = zeros(1, T_len);

for i=1:config.Time_Sim+config.Time
    j= i-floor((i-1)/length(T_ver))*length(T_ver);
  T_traj((i-1)*nT+1: i*nT) = T_ver(j)*ones(nT,1);
end
config.T_traj=T_traj;
config.d_lim = [min(T_ver), max(T_ver)];
