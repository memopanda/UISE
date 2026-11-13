clear all 
close all
%define model constants
plant.al_d1 = 0.544;
plant.al_d2 =  2.65e-7;
plant.al_c1 = 14.28;
plant.al_c2 = plant.al_c1*4.87e-7;
plant.al_c3 = 0.01;
plant.d_c   = 1e-3;  % CO2 outside
plant.al_T1 = 4.57e-4;
plant.al_T2 = 1.38e-5;
plant.d_T   = 10;   % Temperature outside
plant.d_r   = 50;   % solar radiation outside
plant.th_0  = 3.55e-9;
plant.th_1  = 5.11e-6;
plant.th_2  = 2.3e-4;
plant.th_3  = 6.29e-4;
plant.th_4  = 5.2e-5;
plant.dt    = 2*60; % 2min
n = 2;              % #states
nu =2;              % #inputs 
nw = 1;             % #disturbances

%% system dynamics
% x = sdpvar(n,1); 
x = sym('x_%d', [1 2]);
w = sym('w');                   % disturbance
u = sym('u_%d', [1 2]);

f_11 = -plant.al_c1 * (plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*x(2)*x(2)+plant.th_2*x(2)-plant.th_3 ) + x(1)-plant.th_4  ) );
f_w1 = plant.al_c2 * 2^(0.1*x(2)-2.5);  


f_fcn = @(x,w) [-plant.al_c1 * plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*x(2)^2+plant.th_2*x(2)-plant.th_3 ) + x(1)-plant.th_4 )+plant.al_c2*w*2^(0.1*x(2)-2.5) ; 0] ;   
% f_fcnd = @(x,w) (x + plant.dt*f_fcn(x,w) );  % Euler method 
B_fcn = @(x) [plant.al_c3*(plant.d_c -x(1)), 0; -plant.al_T2*(x(2)-plant.d_T),  plant.al_T1 ];
% B_perp = [eye(4); zeros(2,4)];
% Bw = [zeros(1,3),cosx(x(3)),-sinx(x(3)),0]'; 
% Bw_fcn = @(x)[zeros(1,3),cos(x(3)),-sin(x(3)),0]';
% df_11_sub_dx = jacobian(f_11_sub,x);
% df_11_dx = -plant.al_c1 * (plant.th_0*plant.d_r*f_11_sub*(f_11_sub+2*plant.th_0*plant.d_r))/((f_11_sub+plant.th_0*plant.d_r)*(f_11_sub+plant.th_0*plant.d_r)); 
df_11_dx = jacobian(f_11, x);
df_w1_dx = jacobian(f_w1, x);
A = eye(n) + plant.dt*[df_11_dx+ df_w1_dx*w; 0 0] + plant.dt*[ -plant.al_c3*u(1) 0;  0  -plant.al_T2*u(1)];  % I+ T* df/dx
B = plant.dt*[plant.al_c3*(plant.d_c -x(1)), 0; -plant.al_T2*(x(2)-plant.d_T),  plant.al_T1 ];
Bw = plant.dt*[plant.al_c2*2^(0.1*x(2)-2.5)  ;0];
F_fcn_d= @(x,u,w)(x+plant.dt*f_fcn(x,w)+plant.dt*B_fcn(x)*u);



%% for testing with(Chebyshev polynomials) approximated fcns 
% Bw_approx_fcn = @(x)[zeros(1,3),cosx(x(3)),-sinx(x(3)),0]';
% f_vx_approx_fcn = @(x) x(6)*x(5)-plant.g*sinx(x(3));
% Bw_vx_approx_fcn = @(x) cosx(x(3));
% % approximated f_approx, df_dx_fcn
% x = x_store;
% s = sdisplay(f);
% s2 = sdisplay(df_dx);
% s3 = sdisplay(dBw_dx);
% syms x [n 1]
% syms f_approx_fcn [n 1]
% syms df_dx_approx_fcn [n n]
% syms dBw_dx_approx_fcn [n n]
% for i=1:n    
%     f_approx_fcn(i,1) = eval(s{i});    
%     for j=1:n
%         df_dx_approx_fcn(i,j) = eval(s2{i,j});
%         dBw_dx_approx_fcn(i,j) =eval(s3{i,j}); 
%     end
% end
% f_approx_fcn = matlabFunction(f_approx_fcn,'Vars',{x});
% df_dx_approx_fcn = matlabFunction(df_dx_approx_fcn,'Vars',{x});
% dBw_dx_approx_fcn = matlabFunction(dBw_dx_approx_fcn,'Vars',{x});
% x = x_store;
%-------------------------------------------------------------------------


plant.A = matlabFunction(A,'vars', {[x , u, w]} );
plant.B = matlabFunction(B, 'vars', {[x, u, w]});
plant.Bw = matlabFunction(Bw, 'vars',{[x, u, w]});
plant.C = eye(n);
plant.D = eye(nu);
plant.dynamics = F_fcn_d;
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
rccm.L_poly = [x(2)];
% rccm.w_lim  = [ 0.1 0.3];
rccm.w_lim  = [0.1 0.3];
rccm.x_lim  = [0.0005 0.02; 5 40];
% rccm.x_lim  = [0.01 0.02; 20 30];
rccm.x_np   = [5 7];
% rccm.x_np   = [3 3];
rccm.u_lim  = [0 1; 0 1];
rccm.u_np   = [3 3];
rccm.al_lim = [0.005 0.01];
rccm.al_np = 5;


