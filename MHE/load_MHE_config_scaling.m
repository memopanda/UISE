

N_min = max(rccm.N_min, rccm_s.N_min);
config.N_min=ceil(N_min);  % theoritical minimal horizon steps

MHE_d_lim= [-1000   1000 ; rccm.d_lim(2:end,:)];

%%% solver settings
% options = struct;
% options.ipopt.tol = 1e-8;  % default 1e-8
% options.ipopt.dual_inf_tol = 0.001; % 0.01; % default 1;
% options.ipopt.constr_viol_tol =1e-9; % 1e-8; % default 0.0001 
% options.ipopt.compl_inf_tol = 1e-9; %1e-8; % default 0.0001
% options.ipopt.acceptable_constr_viol_tol = 1e-8; %  1e-7;  % default 0.01
% options.ipopt.acceptable_compl_inf_tol =1e-8; % 1e-7;    % default 0.01
% options.ipopt.max_iter = 1200;
% options.ipopt.print_level = 1;  % default 5
% options.print_time = 0; % activate to turn off the NLP output

% options = struct;
% options.ipopt.tol = 1e-8;  % default 1e-8
% options.ipopt.dual_inf_tol = 0.01; % 0.01; % default 1;
% options.ipopt.constr_viol_tol = 1e-5; % 1e-8; % default 0.0001 
% options.ipopt.compl_inf_tol = 1e-4; %1e-8; % default 0.0001
% options.ipopt.acceptable_constr_viol_tol = 0.001; %  1e-7;  % default 0.01
% options.ipopt.acceptable_compl_inf_tol = 0.001; % 1e-7;    % default 0.01
% options.ipopt.max_iter = 2000;
% options.ipopt.print_level = 1;  % default 5
% options.print_time = 0; % activate to turn off the NLP output


% %if tsMHE_estimator
% options = struct;
% options.ipopt.tol = 1e-8;  % default 1e-8
% options.ipopt.acceptable_tol = 1e-6; % default 1e-6
% options.ipopt.acceptable_iter = 25; % default 15
% options.ipopt.dual_inf_tol = 0.000001;  % 0.01; % default 1;
% options.ipopt.constr_viol_tol =1e-9;  % 1e-8; % default 0.0001 
% options.ipopt.compl_inf_tol = 1e-9;   %1e-8; % default 0.0001
% options.ipopt.acceptable_constr_viol_tol = 1e-7;   %1e-7;  % default 0.01
% options.ipopt.acceptable_compl_inf_tol = 1e-7; %1e-7;    % default 0.01
% options.ipopt.max_iter = 2500;
% options.ipopt.print_level = 1;  % default 5
% options.print_time = 0; % activate to turn off the NLP output
% % else



% if  tsMHE_estimator~=2
options = struct;
options.ipopt.tol = 1e-8;  % default 1e-8
options.ipopt.acceptable_tol = 1e-6; % default 1e-6
options.ipopt.acceptable_iter = 15; % default 15
options.ipopt.dual_inf_tol = 0.1; % default 1;
options.ipopt.constr_viol_tol =  1e-5; % default 0.0001 
options.ipopt.compl_inf_tol = 1e-5; % default 0.0001
options.ipopt.acceptable_constr_viol_tol =  1e-3;  % default 0.01
options.ipopt.acceptable_compl_inf_tol =  1e-3;    % default 0.01
options.ipopt.max_iter = 2000;
options.ipopt.print_level = 1;  % default 5
options.print_time = 0; % activate to turn off the NLP output
% 
% end
%% ts-MHE
if tsMHE_estimator==1
config.f_dyn = plant.f_fcn_z_dt;
config.h_dyn = plant.y_fcn_z_dt;

config.T = @(z,y,v)([y(1)-v(1);z(1);z(2)-be*(y(1)-v(1))]); 
config.T_mx = [eye(2) zeros(2,1); -be 0 1]; %  T_mx*[y; z]  = x_est  

config.Txz = @(x,w)([w(1);w(2)+be*x]);  % T_{\sharp} x\to z
 % cost = cost_init + \sum_{i=1}^{Nk}  alpha^{i}  (c_y\|y_{k-i}-y_m_{k-i}\| +c_v \|d_{k-i}\|)
config.cost_init = @(Nk, x_init, x)((rccm_s.alpha_opt^Nk)*2*rccm_s.cx_opt*(x-x_init)'*(x-x_init) );
config.cost_dy = @( y_m, y, d)(rccm_s.cy_opt* (y-y_m(2))'*(y-y_m(2))+ 2*rccm_s.cv_opt*d'*d );  

config.nd =  plant.nv;
config.ny_m = plant.ny;  % # measured output
config.ny = plant.ny-1; % # estimated output
config.nx = rccm_s.n;   % # estimated states in new coordinate
config.nu = plant.nu;
config.alpha = rccm_s.alpha_opt;

config.n_y = config.N*config.ny;  
config.n_y_m = config.N*config.ny_m;   
config.n_d = config.N*config.nd;    
config.n_x =(config.N+1)*config.nx;  


config.xlim  = [rccm.x_lim ; rccm.w_lim] ;
config.xlim_t = rccm.z_lim;
config.ylim  = rccm.y_lim(2,:);
config.dlim  = rccm.d_lim(2:end,:);

% other parameters
config.utraj = rccm.u_lim(:,1);  % kron(ones(1, config.N), rccm.u_lim(:,1));
config.n_u = config.nu;  %config.N*config.nu;
config.n_Nk = 1;
config.n_factor = config.N;


elseif tsMHE_estimator==0
% benchmark MHE

config.f_dyn = @(x,d,y,u)(plant.f_fcn_full_MHE(x,d,u)); 
config.h_dyn = @(x,d,y)(plant.y_fcn_full_MHE(x,d));  

% % the true system only with noises 
%  config.f_dyn = @(x,d,y,u)(plant.fulldyn_d(x,u)); 
%  config.h_dyn = @(x,d,y)(plant.fulloutput_d(x,d(2:3)));  

%%
config.T = @(x,y,v)(x);  %
 
 % cost = cost_init + \sum_{i=1}^{Nk}  alpha^{i}  (c_y\|y_{k-i}-y_m_{k-i}\| +c_v \|d_{k-i}\|)
 
config.cost_init = @(Nk, x_init, x)((rccm.alpha_opt^Nk)*2*rccm.cx_opt*(x-x_init)'*(x-x_init));
config.cost_dy = @( y_m, y, d)( rccm.cy_opt*(y-y_m)'*(y-y_m)+ 2*rccm.cv_opt*(d(2:3)'*d(2:3)+(d(1)-1)'*(d(1)-1) )); 
 
%%
config.nd = plant.nd;  % plant.nv
config.ny_m = plant.ny;  % measured output
config.ny = plant.ny; %plant.ny-1 estimated output
config.nx = rccm.n;    %rccm_s.n
config.nu = plant.nu;
config.alpha = rccm.alpha_opt;

config.n_y = (config.N+1*onoffN)*config.ny;  
config.n_y_m = (config.N+1*onoffN)*config.ny_m;   
config.n_d = (config.N+1*onoffN)*config.nd; 

config.n_x =(config.N+1)*config.nx;  

config.xlim  = [rccm.x_lim ; rccm.w_lim] ;  %  rccm.w_lim 
config.xlim_t = config.xlim;  %   rccm.z_lim
config.ylim  = rccm.y_lim;  % rccm.y_lim(2,:);

config.dlim  =MHE_d_lim;  % rccm.d_lim(2:end,:);



%%
% parameters
config.utraj = rccm.u_lim(:,1);  % kron(ones(1, config.N), rccm.u_lim(:,1));
config.n_u = config.nu;  %config.N*config.nu;
config.n_Nk = 1;
config.n_factor = config.N;
else   %% MHE-based UIE
config.f_dyn = @(x,d,y,u)(plant.f_fcn_full_MHE(x,d,u)); 
config.h_dyn = @(x,d,y)(plant.y_fcn_full_MHE(x,d));  

% % the true system only with noises 
%  config.f_dyn = @(x,d,y,u)(plant.fulldyn_d(x,u)); 
%  config.h_dyn = @(x,d,y)(plant.fulloutput_d(x,d(2:3)));  

%%
config.T = @(x,y,v)(x);  %
 
 % cost = cost_init + \sum_{i=1}^{Nk}  alpha^{i}  (c_y\|y_{k-i}-y_m_{k-i}\| +c_v \|d_{k-i}\|)
 
config.cost_init = @(Nk, x_init, x)((rccm_f.alpha_opt^Nk)*2*rccm_f.cx_opt*(x-x_init)'*(x-x_init));
config.cost_dy = @( y_m, y, d)( rccm_f.cy_opt*(y-y_m)'*(y-y_m)+ 2*rccm_f.cv_opt*(d(2:3)'*d(2:3) )); 
 
%%
config.nd = plant.nd;  % plant.nv
config.ny_m = plant.ny;  % measured output
config.ny = plant.ny; %plant.ny-1 estimated output
config.nx = rccm_f.n;    %rccm_s.n
config.nu = plant.nu;
config.alpha = rccm_f.alpha_opt;

config.n_y = (config.N+1)*config.ny;  
config.n_y_m = (config.N+1)*config.ny_m;   
config.n_d = (config.N+1)*config.nd;    
config.n_x =(config.N+1)*config.nx;  

config.xlim  = [rccm.x_lim ; rccm.w_lim]  ;  %  rccm.w_lim 
config.xlim_t = config.xlim;  %   rccm.z_lim
config.ylim  = rccm.y_lim;  % rccm.y_lim(2,:);

config.dlim  = MHE_d_lim;  % rccm.d_lim(2:end,:);

%%
% parameters
config.utraj = rccm.u_lim(:,1);  % kron(ones(1, config.N), rccm.u_lim(:,1));
config.n_u = config.nu;  %config.N*config.nu;
config.n_Nk = 1;
config.n_factor = config.N;


end


