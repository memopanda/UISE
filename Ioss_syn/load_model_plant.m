%% two types of crops x=[co2, crop1, crop2], u=[tempreture, u_v, u_c]; varying plant.d_T;
%define model constants
% clear all

plant.al_d1 = 0.544;
plant.al_d2 = 2.65e-7;
plant.al_d3 =  0.8;   %0.57;    % 2.type of crop
plant.al_d4 =  1.85e-7;  %    %2.5e-7; % 2.type of crop !!!!!!!!!!!!!!!!
plant.al_p1 = 53;
plant.al_p2 = 45;   % 2.type of crop
plant.al_c1 = 1/4.1;  %14.28;
plant.al_c2 = plant.al_c1*4.87e-7;
plant.al_c3 = plant.al_c1*7.5e-3;
plant.al_c4 = plant.al_c1*1.2e-6;
plant.d_c   = 6e-4;  % CO2 outside
plant.al_T1 = 150/30000;
plant.al_T2 = 1290*7.5e-3/30000;
plant.al_T3 = 0; % 0.2/30000; 
%plant.d_T   = 10;   % Temperature outside  d: known disturbance 
plant.d_r   = 100;   % solar radiation [0 600]  
plant.th_0  = 3.55e-9;
plant.th_1  = 5.11e-6;
plant.th_2  = 2.3e-4;
plant.th_3  = 6.29e-4;
plant.th_4  = 5.2e-5;
plant.dt    = 60; % 1min
plant.beta_k = 1e-8; 
%plant.DT    = 10*plant.dt;  %  min
n = 1;              % # co2
nw = 2;             % # biomasses
nu =3;              % #known inputs  u=[tempreture, u_v, u_c] 
nd = 3;             % #unknown inputs  % [phi_photon (i.e., f_11), measurement/process noises]
ny = 2;             % #outputs
nv = 2;             % # measurement/process noises
% scaling variables: scal*x+shif = x_new (x= (1/scal)(x_new-shif))
% x = S_xv*x_new-b_x; v = (1/scal)*v_new; d = S_d d_new-b_d
scal =1000;
shif = -85;
shif_w = -960;
S_xv = diag(ones(1,n+nw)/scal);  
b_x = [0;  ones(nw,1)*shif/scal];
S_d = diag(ones(1,nd)/scal); 
b_d = [shif_w/scal; 0 ;0];
plant.T_so_xw = @(x)(S_xv*x-b_x);  % de-scaling the states
plant.T_os_xw = @(x)(S_xv\(x+b_x));  % scaling the original states
plant.T_os_d = @(d)(S_d\(d+b_d));
plant.s_x=1/scal;
plant.shif = shif;
plant.shif_w = shif_w;
%% system dynamics for i-ioss (benchmark MHE)
 
x = sym('x_%d', [1 n]);                   % x_1: co2
w = sym('w_%d', [1 nw]);                  % x_2,x_3: two types biomasses
d = sym('d_%d', [1 nd]);                   % unknown inputs (w; v_1; v_2)
u = sym('u_%d', [1 nu]);                   % known inputs


f_fcn_z_2_sym = 2^(0.1*u(1)-2.5);

f_fcn_11_sym = plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*u(1)^2+plant.th_2*u(1)-plant.th_3 ) + x(1)-plant.th_4 ) ;
f_fcn_sym = -plant.al_c1 *plant.dt*( (1-exp(-plant.al_p1*w(1)))*f_fcn_11_sym + d(1))  +plant.al_c2*plant.dt*(w(1)+w(2))*f_fcn_z_2_sym;  % first part of (CO2 ) 
B_fcn_sym = plant.al_c3*plant.dt*(plant.d_c -x(1))*u(2)+  plant.al_c4*plant.dt*u(3); % continous model, input matrix u=[valve;  co2 supplier]
%F_fcn_d= @(w,x,u)(x+plant.dt*f_fcn_sym+plant.dt*B_fcn_sym); % discrete model(Euler method) dyn of CO2 



f_fcn_w_mhe_sym = [plant.dt*plant.al_d1*(1-exp(-plant.al_p1*w(1)))*f_fcn_11_sym-plant.dt*plant.al_d2*w(1)*f_fcn_z_2_sym;...
               plant.dt*plant.al_d3*d(1)-plant.dt*plant.al_d4*w(2)*f_fcn_z_2_sym]; 


%f_fcn_full_MHE= @(x,w,d,u)([x+plant.dt*f_fcn_mhe(w,u,d)+plant.dt*B_fcn_sym; w+plant.dt*f_fcn_w_mhe(w,u,d)]);


% matrices for differential dynamics
df_x1_dx = jacobian(f_fcn_sym, x);
df_x2_dx = jacobian(B_fcn_sym, x);
temp_dx = df_x1_dx+df_x2_dx;
df_x1_dw = jacobian(f_fcn_sym, w);
df_x2_dw = jacobian(B_fcn_sym, w);
temp_dw = df_x1_dw+df_x2_dw;

df_w_dx = jacobian(f_fcn_w_mhe_sym, x);
df_w_dw = jacobian(f_fcn_w_mhe_sym, w);

df_x1_dd = jacobian(f_fcn_sym, d);
df_x2_dd = jacobian(B_fcn_sym, d);
temp_dd = df_x1_dd + df_x2_dd;

df_w_dd = jacobian(f_fcn_w_mhe_sym, d);


temp = [temp_dx, temp_dw; df_w_dx,  df_w_dw];

A = eye(n+nw)+ temp;  % I+ T* [df/dx  df/dw]
B = [temp_dd ; df_w_dd];    % T*[df/dd]
C = [1 0 0; 0 1 1]; % dh/x dh/dw
D = [0 1 0; 0 0 1];  % dh/dd


%% systems dynamics for  i-ioss (two stages MHE)
z = sym('z_%d', [1 nw]);
v = sym('v_%d', [1 nv]);
y = sym('y_1'); %the second output does not appear in the dynamic!
ym = sym('y_2');

be = plant.al_d3/plant.al_c1;



f_fcn_z_11_sym = plant.th_0*plant.d_r*(y(1)-v(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*u(1)^2+plant.th_2*u(1)-plant.th_3 ) + y(1)-v(1)-plant.th_4 ) ;
f_fcn_z_3_sym = plant.dt*plant.al_d1*(1-exp(-plant.al_p1*z(1))) * f_fcn_z_11_sym - plant.dt* plant.al_d2*z(1)*f_fcn_z_2_sym;  % crop 1


f_fcn_z_4_sym = plant.dt*( -plant.al_d3 *(1-exp(-plant.al_p1*z(1))) * f_fcn_z_11_sym+  (be*plant.al_c2*(z(1)+z(2)) -plant.al_d4*z(2))*f_fcn_z_2_sym ); % first part of crop 2
f_fcn_z_5_sym = plant.dt*(plant.al_d4*be-be^2*plant.al_c2) *(y(1)-v(1))*f_fcn_z_2_sym + plant.dt*be*( plant.al_c3*u(2)*(plant.d_c-y(1)+v(1)) + plant.al_c4 * u(3)  );  % second part of crop 2

df_z1_dz= jacobian(f_fcn_z_3_sym, z);
df_z1_dv = jacobian(f_fcn_z_3_sym, v);

df_z2_1_dz = jacobian(f_fcn_z_4_sym,z);
df_z2_2_dz = jacobian(f_fcn_z_5_sym,z);
temp_dz = df_z2_2_dz + df_z2_1_dz;
df_z2_1_dv = jacobian(f_fcn_z_4_sym,v);
df_z2_2_dv = jacobian(f_fcn_z_5_sym,v);
temp_dv = df_z2_2_dv + df_z2_1_dv;


% f_fcn_z_dt = @(z,v,y,u)([z(1)+plant.dt* plant.al_d1*f_fcn_w_3(z,y,v,u)- plant.dt*plant.al_d2*z(1)*f_fcn_w_2(u);...
%                 z(2)+plant.dt*( -plant.al_d3 *f_fcn_w_3(z,y,v,u) + ( be*plant.al_c2*(z(1)+z(2)) -plant.al_d4*z(2)+ (plant.al_d4*be-be^2*plant.al_c2) *(y(1)-v(1)) )*f_fcn_w_2(u) )+plant.dt*be*( plant.al_c3*u(2)*(plant.d_c-y(1)+v(1)) + plant.al_c4 * u(3)  )]);


y_fcn_z_sym =  z(1)+z(2)-be*(y(1)-v(1))+v(2) ;

A_s = eye(nw)+[df_z1_dz; temp_dz];
B_s = [df_z1_dv; temp_dv];
C_s = jacobian(y_fcn_z_sym, z);
D_s = jacobian(y_fcn_z_sym, v);



%% system dynamics for i-IOSS (full order MHE-based UIE)
x = sym('x_%d', [1 n]);                   % x_1: co2
w = sym('w_%d', [1 nw]);                  % x_2,x_3: two types biomasses
d = sym('d_%d', [1 nd]);                   % unknown inputs (w; v_1; v_2)
u = sym('u_%d', [1 nu]);                   % known inputs
vf = sym('uf_%d', [1 nd-1]);                 % future unknown inputs (v_1,v_2)



f_fcn_z_2_sym = 2^(0.1*u(1)-2.5);

f_fcn_11_sym = plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*u(1)^2+plant.th_2*u(1)-plant.th_3 ) + x(1)-plant.th_4 ) ;
f_fcn_sym = -plant.al_c1 *plant.dt*( (1-exp(-plant.al_p1*w(1)))*f_fcn_11_sym + d(1))  +plant.al_c2*plant.dt*(w(1)+w(2))*f_fcn_z_2_sym;  % first part of (CO2 ) 
B_fcn_sym = plant.al_c3*plant.dt*(plant.d_c -x(1))*u(2)+  plant.al_c4*plant.dt*u(3); % continous model, input matrix u=[valve;  co2 supplier]
%F_fcn_d= @(w,x,u)(x+plant.dt*f_fcn_sym+plant.dt*B_fcn_sym); % discrete model(Euler method) dyn of CO2 



f_fcn_w_mhe_sym = [plant.dt*plant.al_d1*(1-exp(-plant.al_p1*w(1)))*f_fcn_11_sym-plant.dt*plant.al_d2*w(1)*f_fcn_z_2_sym;...
               plant.dt*plant.al_d3*d(1)-plant.dt*plant.al_d4*w(2)*f_fcn_z_2_sym]; 


%f_fcn_full_MHE= @(x,w,d,u)([x+plant.dt*f_fcn_mhe(w,u,d)+plant.dt*B_fcn_sym; w+plant.dt*f_fcn_w_mhe(w,u,d)]);


% matrices for differential dynamics
df_x1_dx = jacobian(f_fcn_sym, x);
df_x2_dx = jacobian(B_fcn_sym, x);
temp_dx = df_x1_dx+df_x2_dx;
df_x1_dw = jacobian(f_fcn_sym, w);
df_x2_dw = jacobian(B_fcn_sym, w);
temp_dw = df_x1_dw+df_x2_dw;

df_w_dx = jacobian(f_fcn_w_mhe_sym, x);
df_w_dw = jacobian(f_fcn_w_mhe_sym, w);

df_x1_dd = jacobian(f_fcn_sym, d);
df_x2_dd = jacobian(B_fcn_sym, d);
temp_dd = df_x1_dd + df_x2_dd;

df_w_dd = jacobian(f_fcn_w_mhe_sym, d);


temp = [temp_dx, temp_dw; df_w_dx,  df_w_dw];

A_f = eye(n+nw)+ temp;  % I+ T* [df/dx  df/dw]
B_f = [temp_dd(:,2:end), zeros(n, nd-1) ; df_w_dd(:,2:end), zeros(nw, nd-1)];    % T*[df/dv df/dvf] 
C_0 = [1 0 0; 0 1 1]; % dh/dx dh/dw
C_1 = [A_f(1,:); A_f(2,:)+A_f(1+nw,:)]; % dhf/dx dhf/dw
C_f = [C_0; C_1];   % 
D_0 = [1 0;  0 1];  % [dh/dv ]
D_1 = [D_0 zeros(size(D_0,1), nd-1 ) ];  % [dh/dv dh/dvf]
D_2 = [temp_dd(:,2:end); df_w_dd(1, 2:end)+df_w_dd(2, 2:end)]; % dhf/dv
D_f = [D_1; D_2  eye(nd-1)]; % [dh/dv dh/dvf; dhf/dv dhf/dvf]

%% models for simulations (true full system)

f_fcn_11 = @ (x,u) (plant.th_0*plant.d_r*(x(1)-plant.th_4)  /( plant.th_0*plant.d_r/(-plant.th_1*u(1)^2+plant.th_2*u(1)-plant.th_3 ) + x(1)-plant.th_4 )   );
f_fcn = @(w,x,u,k) (-plant.al_c1 *( 2-exp(-plant.al_p1*w(1)) -exp(-plant.al_p2*w(2))+plant.beta_k*k )*f_fcn_11(x,u) +plant.al_c2*(w(1)+w(2))*2^(0.1*u(1)-2.5)  ) ;  % first part of (CO2 ) 
B_fcn = @(x) [plant.al_c3*(plant.d_c -x(1)),  plant.al_c4]; % continous model, input matrix u=[valve;  co2 supplier]
F_fcn_d= @(w,x,u,k)(x+plant.dt*f_fcn(w,x,u,k)+plant.dt*B_fcn(x)*[u(2);u(3)]); % discrete model(Euler method) dyn of CO2 


f_fcn_w = @(w,x,u)[(1-exp(-plant.al_p1*w(1)))*plant.al_d1* f_fcn_11(x,u)-plant.al_d2*w(1)*2^(0.1*u(1)-2.5);...
          (1-exp(-plant.al_p2*w(2)))*plant.al_d3*f_fcn_11(x,u)-plant.al_d4*w(2)*2^(0.1*u(1)-2.5) ]; % correct dyn of crops
F_fcn_full_d = @(x,u,k)([ F_fcn_d(x(n+1:nw+n),x(1),u,k); x(n+1:n+nw) + plant.dt*f_fcn_w(x(n+1:nw+n),x(1),u)]); % 
F_fcn_full_d_scal = @(x,u,k)(S_xv\F_fcn_full_d( S_xv*x - b_x, u,k)+S_xv\b_x );
y_fcn_full_d = @(x,v)([x(1)+v(1); x(2)+x(3)+v(2)]);


%% model for two stages MHE
% w_1 = z_2 ;
% w_2 = z_3 - be*(y_1-d(2)), x = z_1 = y(1)-d(2);

 
f_fcn_w_2 = @(u)(2^(0.1*u(1)-2.5));
f_fcn_w_3 = @(z,y,v,u)(  (1-exp(-plant.al_p1*z(1)) ) * f_fcn_11(y-v,u) );
f_fcn_w_4 = @(z,y,v)(plant.al_c2*be*(z(1)+z(2)) -plant.al_d4*z(2)+ be*(plant.al_d4-be*plant.al_c2) *(y(1)-v(1)) ); 

f_fcn_z_dt = @(z,v,y,u)([z(1)+plant.dt* plant.al_d1*f_fcn_w_3(z,y,v,u)- plant.dt*plant.al_d2*z(1)*f_fcn_w_2(u);...
                z(2)+plant.dt*( -plant.al_d3 *f_fcn_w_3(z,y,v,u) ) + plant.dt*f_fcn_w_4(z,y,v) *f_fcn_w_2(u)+plant.dt*be*( plant.al_c3*u(2)*(plant.d_c-y(1)+v(1)) + plant.al_c4 * u(3)  )]);
f_fcn_z_dt_scal = @(z,v,y,u)( S_xv(2:3, 2:3)\f_fcn_z_dt( S_xv(n+1:end,n+1:end)*z - b_x(n+1:end), S_xv(1:nv,1:nv)*v, S_xv(1:ny,1:ny)*y,  u) + S_xv(2:3, 2:3)\b_x(n+1:end)    );  % fct only depends on y(1), which is not shifted. 
y_fcn_z_dt = @(z,v,y)(  z(1)+z(2)-be*(y(1)-v(1))+v(2)  );





%% model for benchmark MHE and full-order MHE-based UIE
f_fcn_mhe = @(x, w, u, d) (-plant.al_c1 * ( (1-exp(-plant.al_p1*w(1)))* f_fcn_11(x,u)+ d(1) ) +plant.al_c2*(w(1)+w(2))*f_fcn_w_2(u)  );
f_fcn_w_mhe = @(x, w, u, d)[plant.al_d1*(1-exp(-plant.al_p1*w(1)))*f_fcn_11(x,u)-plant.al_d2*w(1)*f_fcn_w_2(u);...
          plant.al_d3*d(1)-plant.al_d4*w(2)*f_fcn_w_2(u) ]; 
f_fcn_full_MHE= @(x,d,u)([x(1)+plant.dt*f_fcn_mhe(x(1),x(n+1:nw+n),u,d)+plant.dt*B_fcn(x(1))*[u(2);u(3)]; x(n+1:nw+n)+plant.dt*f_fcn_w_mhe(x,x(n+1:nw+n),u,d)]);
f_fcn_full_MHE_scal = @(x,d,u)( S_xv\f_fcn_full_MHE( S_xv*x - b_x, S_d*d-b_d   , u) +S_xv\b_x );
y_fcn_full_MHE = @(x,d)([x(1)+d(2); x(2)+x(3)+d(3)]);


%% model for MHE-based UIE



%%  Settings and Preprocessing

plant.fulldyn_d = F_fcn_full_d_scal;  % true full system (process)
plant.fulloutput_d = y_fcn_full_d;  % true full system (output)

%plant.wk_d = @(x)(  1-exp(-plant.al_p2*(1/scal)*(x(3)-shif)));  % unscaled wk
plant.wk_d = @(x,u,k)( (1-exp(-plant.al_p2*x(3)))*f_fcn_11(x,u)+plant.beta_k*k);  % unscaled wk

plant.f_fcn_z_dt = f_fcn_z_dt_scal;   % model for two stages MHE (process)
plant.y_fcn_z_dt = y_fcn_z_dt;   % model for two stages MHE (output)


plant.f_fcn_full_MHE = f_fcn_full_MHE_scal;   % model for benchmark MHE (process)
plant.y_fcn_full_MHE = y_fcn_full_MHE;     %  model for benchmark MHE (output)




plant.n = 1;              % # co2
plant.nw = 2;             % # biomasses
plant.nu =3;              % #known inputs  u=[tempreture, u_v, u_c] 
plant.nd = 3;             % #unknown inputs  % [phi_photon (i.e., f_11), measurement/process noises]
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
rccm.w_lim  = [ 0.08   0.1; 0.08   0.1];  % [0.001  0.01  0.05   0.1  0.3];
%rccm.w_lim  = [ 0.08   0.083; 0.08   0.083];  % [0.001  0.01  0.05   0.1  0.3];
rccm.w_np   = [4, 4];
rccm.x_lim  = [0 0.00275];
% rccm.x_lim  = [0 0.0015];
rccm.x_np   = [5];

rccm.u_lim  = [25 25; 0 0 ; 0 0 ];  % inner temperature, valve for CO2, supplier CO2
rccm.u_np   = [1 1 1];

%rccm.d_lim = [0.95  1  ;  -5e-6  5e-6; -5e-6 5e-6 ]; % disturbance,meas.noise 1,  meas.noise 2 --paper
rccm.d_lim = [-1000  1000  ;  -5e-6  5e-6; -5e-6 5e-6 ]; % disturbance, meas.noise 1,  meas.noise 2
rccm.d_np  = [2, 2, 2];

rccm.v_lim = rccm.d_lim(2:3,:);
rccm.v_np = rccm.d_np(2:3);

rccm.y_lim = [rccm.x_lim+rccm.d_lim(2,:); rccm.w_lim(1,:)*2+rccm.d_lim(3,:)];
rccm.y_np  = [4 4];

rccm.z_lim = [rccm.w_lim(1,:); rccm.w_lim(2,:)+be*rccm.x_lim];
rccm.z_np  = [rccm.w_np(1),2];



% configurations for ioss synthesis
rccm.al_lim = [0.36 0.36];  %0.48
rccm.al_np = 1;

rccm.P_eps = 1e-5;  %1e-5
rccm.cvy_eps = 1e9;  % 1e9
% rccm.options=[1e-4 400 0 15 1];
rccm.options=[1e-6 500 0 20 1];


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



