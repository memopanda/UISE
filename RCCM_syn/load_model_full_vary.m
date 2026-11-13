%% only one type of crop; varying plant.d_T
%define model constants
plant.al_d1 = 0.544;
plant.al_d2 =  2.65e-7;
plant.al_c1 = 14.28;  %14.28;
plant.al_c2 = plant.al_c1*4.87e-7;
plant.al_c3 = 0.01;
plant.al_c4 = 1.8e-5;   % co2 supplier
plant.d_c   = 7e-4;  % CO2 outside
plant.al_T1 = 4.57e-4;
plant.al_T2 = 1.38e-5;
% plant.d_T   = 10;   % 10  Temperature outside
plant.d_r   = 100;   % 100 solar radiation outside

plant.th_0  = 3.55e-9;
plant.th_1  = 5.11e-6;
plant.th_2  = 2.3e-4;
plant.th_3  = 6.29e-4;
plant.th_4  = 5.2e-5;
plant.dt    = 60; % 2min
plant.DT    = 15*plant.dt;  % 30min 
n = 2;              % #states of fast submodel
nu =3;              % #inputs 
nw = 1;             % #disturbances (#states of other slow model)

%% system dynamics for RCCM
% x = sdpvar(n,1); 
x = sym('x_%d', [1 n]);
w = sym('w');                   % disturbance
u = sym('u_%d', [1 nu]);
d = sym('d');

f_11 = -plant.al_c1 * (plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*x(2)*x(2)+plant.th_2*x(2)-plant.th_3 ) + x(1)-plant.th_4  ) );  %phi_phot ohne w % continous model
f_w1 = plant.al_c2 * 2^(0.1*x(2)-2.5);  % phi_resp ohne w % continous model
f_w2 = 1-exp(-53*w); % continous model

f_fcn = @(w,x) [-plant.al_c1 * plant.th_0*plant.d_r*(1-exp(-53*w))*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*x(2)^2+plant.th_2*x(2)-plant.th_3 ) + x(1)-plant.th_4 )+plant.al_c2*w*2^(0.1*x(2)-2.5); 0] ;  % continous model 
B_fcn = @(x,d) [plant.al_c3*(plant.d_c -x(1)), 0  plant.al_c4; -plant.al_T2*(x(2)-d),  plant.al_T1  0]; % continous model
F_fcn_d= @(w,x,u,d)(x+plant.dt*f_fcn(w,x)+plant.dt*B_fcn(x,d)*u); % discrete fast submodel(Euler method)
% matrices for differential dynamics
df_11_dx = jacobian(f_11, x);
df_w1_dx = jacobian(f_w1, x);
df_w2_dw = jacobian(f_w2, w);

A = eye(n) + plant.dt*[df_11_dx*f_w2+ df_w1_dx*w; 0 0] + plant.dt*[ -plant.al_c3*u(1) 0;  0  -plant.al_T2*u(1)];  % I+ T* df/dx
B = plant.dt*[plant.al_c3*(plant.d_c -x(1)), 0  plant.al_c4; -plant.al_T2*(x(2)-d),  plant.al_T1  0 ];
Bw = plant.dt*[f_w1+f_11*df_w2_dw  ;0];


%% Full system 
plant.n_f = n+nw;  % #states of full system;

f_fcn_w11= @(x)(plant.al_d1 * plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*x(2)^2+plant.th_2*x(2)-plant.th_3 ) + x(1)-plant.th_4 )  );
f_fcn_w = @(w,x)((1-exp(-53*w))*f_fcn_w11(x)-plant.al_d2*w*2^(0.1*x(2)-2.5));  
F_fcn_full_d = @(w, x,u,d)([w;x]+ plant.dt*[f_fcn_w(w,x); f_fcn(w,x)]+plant.dt*[zeros(1,nu); B_fcn(x,d)]*u); 
f_fcn_w_dt = @(w,x)(w+plant.dt*f_fcn_w(w,x));  % discrete slow submodel (orginal)
f_fcn_w_DT = @(w,x)(w+plant.DT*f_fcn_w(w,x));   % discrete slow submodel (coarse)
H_fcn_w_DT = @(x)(plant.DT*[f_fcn_w11(x) -plant.al_d2*2^(0.1*x(2)-2.5)]);
g_sym_w_DT = [1-exp(-53*w); w ];
dg_sym_w_DT = jacobian(g_sym_w_DT, w);
% % debug
% w0=0.2; x0=[0.01; 20]; 
% tt1= f_fcn_w_DT(w0, x0 )-w0;
% tt2 =H_fcn_w_DT(x0)*g_fcn_w_DT(w0);
% if abs(tt1-tt2)<1e-8 
%     disp('no error!')
% end 

%%

plant.fulldyn_d = F_fcn_full_d;  % original full model
plant.f_fcn_w_dt = f_fcn_w_dt;   % original slow submodel
plant.f_fcn_w_DT = f_fcn_w_DT;     % coarse slow submodel
plant.H_fcn_w_DT = H_fcn_w_DT;
plant.dg_fcn_w_DT = matlabFunction(dg_sym_w_DT);

plant.A = matlabFunction(A,'vars', {[x , u, w, d]} );
plant.B = matlabFunction(B, 'vars', {[x, u, w, d]});
plant.Bw = matlabFunction(Bw, 'vars',{[x, u, w, d]});
plant.C = eye(n);
plant.D = eye(nu);
plant.dynamics = F_fcn_d; % original fast submodel
plant.w = w;

plant.x = x;
plant.u = u;

plant.n = n; plant.nu=nu; plant.nw = nw; 


% if controller.opt_pos_dev == 0
%     C = [eye(n); zeros(nu,n)]; D =  [zeros(n,nu); 1*eye(nu)];  % consider both states and inputs
% elseif controller.opt_pos_dev == 1
%     C = [eye(2) zeros(2,4);zeros(nu,n)]; 
%     D = [zeros(2,nu); 1*eye(nu)];                               % consider only position states and inputs
% end
% nz = size(C,1); 
% plant.C = C; plant.D = D; plant.nz = nz;


% %% Optimization vars
% %----------------- settings for searching CCM ----------------------
% x = sdpvar(n,1); x_store = x;
% w = sdpvar(nw,1); % disturbance and uncert param w = [d; theta]. theta is 
% % shifted by the nominal parameter \bar{\theta}_0 = 1/plant.m
% u = sdpvar(nu,1); % need input as variable since G(u) used for synthesis
% 
% W_states_index = [3 4]; %specify which states are used to parametrized the CCM
% 
% % ------------------ Stateconstraints for metric synthesis ----------------
% 
% %specify constraints of as ellipsoidals
% state_set.box_lim = [state_set.p_lim^2-x(3)^2; state_set.vx_lim^2-x(4)^2; ...
%     state_set.pd_lim^2-x(6)^2;  state_set.vz_lim^2-x(5)^2]*0.001; % scaling 0.001???
% %state_set.box_lim = [state_set.box_lim; state_set.d_lim^2-w^2]; 
% state_set.num_consts_4_W_states = 2;        % # constraints from box_lim that involve states on which the metric W depends
% state_set.other_lim_states = [x(6);x(5)]; 
% state_set.lagrange_deg_W = 4;               % for the bounds of W
% state_set.lagrange_deg_ccm = 4;             % for ccm condition; had to reduce since there are 3 new variables: theta, u(1) and u(2)
% 
% u_av = mean([state_set.u_min, state_set.u_max]);
% u_lim = state_set.u_max - u_av;
% state_set.box_lim = [state_set.box_lim; u_lim^2 - (u(1) - u_av)^2; u_lim^2 - (u(2) - u_av)^2]; % add scaling? 

rccm.W_poly = [];  % {'x_2', 'x_2*x_1'}  x(2)*x(2)
rccm.L_poly = [x(1) x(2)  ];

rccm.w_lim  = [0.1  0.3];  %  1-exp(-53*w) almost =1
rccm.x_lim  = [0.0001 0.018; 5 40];
rccm.x_np   = [6 7];
% rccm.x_np   = [3 3];
rccm.u_lim  = [0 1; 0 1; 0 1];
rccm.u_np   = [2 2 2];

rccm.d_lim = [5 10 15];

rccm.al_lim = [0.4 0.5];
rccm.al_np = 2;
rccm.ccm_eps = 1e-2;
rccm.c_x1=10; 
rccm.c_u=0.0005;
rccm.options=[1e-4 300 0 15 1];

