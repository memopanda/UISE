%% only one type of crop;  constant plant.d_T
%define model constants
plant.al_c1 = 0.2;
plant.al_c2 =  0.5;
plant.al_c3 = 0.2;  %14.28;
plant.al_c4 = 0.1;
plant.al_d1 = 2;

plant.con_k  = 0.1; %0.1;
plant.dt    = 0.25; %

n = 1;              % # x_1
nw = 2;             % # x_2, x_3
nu =1;              % #known inputs   
nd = 3;             % #unknown inputs  % [phi_photon (i.e., f_11), measurement/process noises]
ny = 2;             % #outputs
nv = 2;             % # measurement/process noises

%% system dynamics for i-ioss (benchmark MHE)
 
x = sym('x_%d', [1 n]);                   % x_1: 
w = sym('w_%d', [1 nw]);                  % x_2,x_3: 
d = sym('d_%d', [1 nd]);                   % unknown inputs (d_1; v_1; v_2)
u = sym('u_%d', [1 nu]);                   % known inputs



f_fcn_sym = plant.dt*(plant.al_c1 *w(1)-plant.al_c2*x(1)*w(2)-plant.al_d1*d(1)+plant.al_d1*plant.al_c4*w(2) );  

f_fcn_w_mhe_sym = [plant.dt*(plant.al_c2*x(1)*w(2)-plant.al_c1*w(1)); plant.dt*(plant.al_c1*w(1)-plant.al_c2*x(1)*w(2)+d(1)-plant.al_c4*w(2))  ];


% matrices for differential dynamics
temp_dx = jacobian(f_fcn_sym, x);
temp_dw =jacobian(f_fcn_sym, w);

df_w_dx = jacobian(f_fcn_w_mhe_sym, x);
df_w_dw = jacobian(f_fcn_w_mhe_sym, w);

df_x1_dd = jacobian(f_fcn_sym, d);


df_w_dd = jacobian(f_fcn_w_mhe_sym, d);


temp = [temp_dx, temp_dw; df_w_dx,  df_w_dw];

A = eye(n+nw)+ temp;  % I+ T* [df/dx  df/dw]
B = [df_x1_dd ; df_w_dd];    % T*[df/dd]
C = [1 0 0; 1 1 1]; % dh/x dh/dw
D = [0 1 0; 0 0 1];  % dh/dd



%% systems dynamics for  i-ioss (two stages MHE)
z = sym('z_%d', [1 nw]);
v = sym('v_%d', [1 nv]);
y = sym('y_1'); %the second output does not appear in the dynamic!
ym = sym('y_2');

be = 1/plant.al_d1;


f_fcn_z_temp = plant.al_c2*(y(1)-v(1))*(z(2)-be*y(1)+be*v(1));
f_fcn_z_3_sym = plant.dt*(f_fcn_z_temp-plant.al_c1*z(1)); 


f_fcn_z_4_sym = -(1+be)*plant.dt*( f_fcn_z_temp-plant.al_c1*z(1) ); 


df_z1_dz= jacobian(f_fcn_z_3_sym, z);
df_z1_dv = jacobian(f_fcn_z_3_sym, v);

temp_dz = jacobian(f_fcn_z_4_sym,z);

temp_dv = jacobian(f_fcn_z_4_sym,v);



% f_fcn_z_dt = @(z,v,y,u)([z(1)+plant.dt* plant.al_d1*f_fcn_w_3(z,y,v,u)- plant.dt*plant.al_d2*z(1)*f_fcn_w_2(u);...
%                 z(2)+plant.dt*( -plant.al_d3 *f_fcn_w_3(z,y,v,u) + ( be*plant.al_c2*(z(1)+z(2)) -plant.al_d4*z(2)+ (plant.al_d4*be-be^2*plant.al_c2) *(y(1)-v(1)) )*f_fcn_w_2(u) )+plant.dt*be*( plant.al_c3*u(2)*(plant.d_c-y(1)+v(1)) + plant.al_c4 * u(3)  )]);


y_fcn_z_sym =  z(1)+z(2)+be*(y(1)-v(1))+v(2) ;

A_s = eye(nw)+[df_z1_dz; temp_dz];
B_s = [df_z1_dv; temp_dv];
C_s = jacobian(y_fcn_z_sym, z);
D_s = jacobian(y_fcn_z_sym, v);


%% system dynamics for i-IOSS (full order MHE-based UIE)
x = sym('x_%d', [1 n]);                   % x_1: 
w = sym('w_%d', [1 nw]);                  % x_2,x_3: 
d = sym('d_%d', [1 nd]);                   % unknown inputs (d_1; v_1; v_2)
u = sym('u_%d', [1 nu]);                   % known inputs
vf = sym('uf_%d', [1 nd-1]);                 % future unknown inputs (v_1,v_2)

f_fcn_sym = plant.dt*(plant.al_c1 *w(1)-plant.al_c2*x(1)*w(2)-plant.al_d1*d(1)+plant.al_d1*plant.al_c4*w(2) );  

f_fcn_w_mhe_sym = [plant.dt*(plant.al_c2*x(1)*w(2)-plant.al_c1*w(1)); plant.dt*(plant.al_c1*w(1)-plant.al_c2*x(1)*w(2)+d(1)-plant.al_c4*w(2))  ];



% matrices for differential dynamics
temp_dx = jacobian(f_fcn_sym, x);
temp_dw =jacobian(f_fcn_sym, w);

df_w_dx = jacobian(f_fcn_w_mhe_sym, x);
df_w_dw = jacobian(f_fcn_w_mhe_sym, w);

df_x1_dd = jacobian(f_fcn_sym, d);


df_w_dd = jacobian(f_fcn_w_mhe_sym, d);


temp = [temp_dx, temp_dw; df_w_dx,  df_w_dw];

% A = eye(n+nw)+ temp;  % I+ T* [df/dx  df/dw]
% B = [df_x1_dd ; df_w_dd];    % T*[df/dd]
% C = [1 0 0; 1 1 1]; % dh/x dh/dw
% D = [0 1 0; 0 0 1];  % dh/dd

A_f = eye(n+nw)+ temp;  % I+ T* [df/dx  df/dw]
B_f = [df_x1_dd(:,2:end), zeros(n, nd-1) ; df_w_dd(:,2:end), zeros(nw, nd-1)];    % T*[df/dv df/dvf] 
C_0 = [1 0 0; 1 1 1]; % dh/dx dh/dw
C_1 = [A_f(1,:); A_f(1,:) + A_f(2,:)+A_f(1+nw,:)]; % dhf/dx dhf/dw
C_f = [C_0; C_1];   % 
D_0 = [1 0;  0 1];  % [dh/dv ]
D_1 = [D_0 zeros(size(D_0,1), nd-1 ) ];  % [dh/dv dh/dvf]
D_2 = [df_x1_dd(:,2:end); df_x1_dd(:, 2:end)+ df_w_dd(1, 2:end)+df_w_dd(2, 2:end)]; % dhf/dv
D_f = [D_1; D_2  eye(nd-1)]; % [dh/dv dh/dvf; dhf/dv dhf/dvf]

%% models for simulations (true full system)


f_fcn =@(w,x,u,k)(plant.al_c1 *w(1)-plant.al_c2*x(1)*w(2)-plant.al_d1*(plant.con_k*k+plant.al_c3)*x(1)^2+plant.al_d1*plant.al_c4*w(2));  

F_fcn_d= @(w,x,u,k)(x+plant.dt*f_fcn(w,x,u,k)); % 


f_fcn_w = @(w,x,u,k)[(plant.al_c2*x(1)*w(2)-plant.al_c1*w(1)); (plant.al_c1*w(1)-plant.al_c2*x(1)*w(2)+(plant.con_k*k+plant.al_c3)*x(1)^2-plant.al_c4*w(2))]; % correct dyn of crops
F_fcn_full_d_scal = @(x,u,k)([ F_fcn_d(x(n+1:nw+n),x(1),u,k); x(n+1:n+nw) + plant.dt*f_fcn_w(x(n+1:nw+n),x(1),u,k)]); % 

y_fcn_full_d = @(x,v)([x(1)+v(1); x(1)+x(2)+x(3)+v(2)]);


%% model for two stages MHE
% w_1 = z_2 ;
% w_2 = z_3 - be*(y_1-d(2)), x = z_1 = y(1)-d(2);

 



f_fcn_z_1 = @(z,y,v) (plant.al_c2*(y(1)-v(1))*(z(2)-be*y(1)+be*v(1)));
f_fcn_z_3 = @(z,y,v) (plant.dt*(f_fcn_z_1(z,y,v)-plant.al_c1*z(1))); 


f_fcn_z_4 = @(z,y,v) (-(1+be)*plant.dt*( f_fcn_z_1(z,y,v)-plant.al_c1*z(1) )); 



f_fcn_z_dt_scal = @(z,v,y,u)([ z(1)+f_fcn_z_3(z,y,v);z(2)+f_fcn_z_4(z,y,v)]  );  % fct only depends on y(1), which is not shifted. 
y_fcn_z_dt = @(z,v,y)(  z(1)+z(2)+be*(y(1)-v(1))+v(2)  );





%% model for benchmark MHE and full-order MHE-based UIE




f_fcn_mhe = @(x, w, u, d) (plant.al_c1 *w(1)-plant.al_c2*x(1)*w(2)-plant.al_d1*d(1)+plant.al_d1*plant.al_c4*w(2)  );
f_fcn_w_mhe = @(x, w, u, d)[plant.al_c2*x(1)*w(2)-plant.al_c1*w(1);...
          plant.al_c1*w(1)-plant.al_c2*x(1)*w(2)+d(1)-plant.al_c4*w(2)]; 
f_fcn_full_MHE_scal= @(x,d,u)([x(1)+plant.dt*f_fcn_mhe(x(1),x(n+1:nw+n),u,d); x(n+1:nw+n)+plant.dt*f_fcn_w_mhe(x,x(n+1:nw+n),u,d)]);

y_fcn_full_MHE = @(x,d)([x(1)+d(2); x(1)+x(2)+x(3)+d(3)]);






%%  Settings and Preprocessing

plant.fulldyn_d = F_fcn_full_d_scal;  % true full system (process)
plant.fulloutput_d = y_fcn_full_d;  % true full system (output)

%plant.wk_d = @(x)(  1-exp(-plant.al_p2*(1/scal)*(x(3)-shif)));  % unscaled wk
plant.wk_d = @(x,u,k)( (plant.con_k*k+plant.al_c3)*x(1)^2 );  % unscaled wk

plant.f_fcn_z_dt = f_fcn_z_dt_scal;   % model for two stages MHE (process)
plant.y_fcn_z_dt = y_fcn_z_dt;   % model for two stages MHE (output)


plant.f_fcn_full_MHE = f_fcn_full_MHE_scal;   % model for benchmark MHE (process)
plant.y_fcn_full_MHE = y_fcn_full_MHE;     %  model for benchmark MHE (output)




plant.n = 1;              % # co2
plant.nw = 2;             % # biomasses
plant.nu =1;              % #known inputs  u=[tempreture, u_v, u_c] 
plant.nd = 3;             % #unknown inputs  % [d, measurement/process noises]
plant.ny = 2;             % #outputs
plant.nv = 2;             % # measurement/process noises

plant.w = w;
plant.d = d;
plant.u = u;
plant.x = x;

plant.z = z;
plant.v = v;
plant.y = [y,ym];



% limiation of variables
rccm.w_lim  = [ 0   4; 0   4];  % [0.001  0.01  0.05   0.1  0.3];
%rccm.w_lim  = [ 0.08   0.083; 0.08   0.083];  % [0.001  0.01  0.05   0.1  0.3];
rccm.w_np   = [2, 4];
rccm.x_lim  = [0 4];
% rccm.x_lim  = [0 0.0015];
rccm.x_np   = [4];

rccm.u_lim  = [0 0 ];  % inner temperature, valve for CO2, supplier CO2
rccm.u_np   = [2 ];

%rccm.d_lim = [0.95  1  ;  -5e-6  5e-6; -5e-6 5e-6 ]; % disturbance,meas.noise 1,  meas.noise 2 --paper
rccm.d_lim = [-1000  1000  ;  -1e-2  1e-2; -1e-2 1e-2 ]; % disturbance, meas.noise 1,  meas.noise 2
rccm.d_np  = [2, 2, 2];

rccm.v_lim = rccm.d_lim(2:3,:);
rccm.v_np = rccm.d_np(2:3);

rccm.y_lim = [rccm.x_lim+rccm.d_lim(2,:); rccm.w_lim(1,:)*3+rccm.d_lim(3,:)];
rccm.y_np  = [4 4];

rccm.z_lim = [rccm.w_lim(1,:); rccm.w_lim(2,:)+be*rccm.x_lim];
rccm.z_np  = [rccm.w_np(1),2];



% configurations for ioss synthesis
rccm.al_lim = [0.8 0.8];  %0.48
rccm.al_np = 1;

rccm.P_eps = 1e-1;  %1e-5
rccm.cvy_eps = 1e2;  % 1e9
 rccm.options=[1e-4 400 0 15 1];
%rccm.options=[1e-6 500 0 20 1];


% i-ioss for benchmark MHE
rccm.A = matlabFunction(A,'vars', {[x, w, d, u]} );
rccm.B = matlabFunction(B, 'vars', {[x, w, d, u]});
rccm.C = C;
rccm.D = D;
rccm.n = n+nw; 
rccm.nv = plant.nd;
rccm.ny = plant.ny;

% i-ioss for MHE-based UIE
rccm_f=rccm;
rccm_f.A = matlabFunction(A_f,'vars', {[x, w, d, u]} );
rccm_f.B = matlabFunction(B_f, 'vars', {[x, w, d, u]});
rccm_f.C = matlabFunction(C_f, 'vars', {[x, w, d, u]});
rccm_f.D = matlabFunction(D_f, 'vars', {[x, w, d, u]});
rccm_f.n = n+nw; 
rccm_f.nv = plant.nd-1+plant.nd-1;
rccm_f.ny = plant.ny*2;


% i-ioss for two stages MHE
rccm_s=rccm;
rccm_s.A = matlabFunction(A_s, 'vars', {[z, v, y, u]});
rccm_s.B = matlabFunction(B_s, 'vars', {[z, v, y, u]});
rccm_s.C = C_s;
rccm_s.D = D_s;
rccm_s.n = nw; 
rccm_s.nv = plant.nv;
rccm_s.ny = plant.ny-1;

%% generate gridding points for MHE and two stages MHE
Gen_gridding_points



