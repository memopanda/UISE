%%  Collect measurements of true system




xw_init = [0.05; 0.5; 0];  % real initial state
xw_pri = [0.5; 1;0.1];  % priori
u_const = [0];  % not used



real_traj = zeros(plant.n+plant.nw,config.N_sim+1);
ym_store = zeros(plant.ny,config.N_sim);
%% debug
v_store=2*v_noise*(rand(plant.nv,config.N_sim)-0.5*ones(plant.nv,config.N_sim)); % --> paper

wk_store = zeros(1, config.N_sim);  % true disturbance

%% Generate measurments and states of true model
real_traj(:,1)=xw_init;
for k= 1: config.N_sim
[xw_new]=plant.fulldyn_d(xw_init(1:plant.nw+plant.n), u_const,k);
real_traj(:,k+1) = xw_new;
ym_store(:,k)=plant.fulloutput_d(xw_init, v_store(:,k));

wk_store(:,k)=plant.wk_d(xw_init,u_const,k);
xw_init = xw_new;

end
%%


% xw_pri_t = xw_init_t;
 d_est = [0; zeros(plant.nv,1)]; 

config.x_init = xw_init;
config.x_init_est = xw_pri; 
config.u_true = u_const;
config.v_true = v_store;
config.w_true = wk_store;
%% MHE Simulation 


% lower and upper bound on constraints
c_lb = [zeros(po.n_eq,1); -inf(po.n_ineq,1)];

c_ub = zeros(po.n_eq + po.n_ineq,1);

% lower and upper bound 
xw_lb = -inf(po.n_var,1);
xw_ub = inf(po.n_var,1);
%  setting bounds for each variable
% xw_lb(end-config.n_w+1:end) = zeros(config.n_w,1);
% xw_ub(end-config.n_e+1:end) = diag([0.001, 0.01])*ones(config.n_e,1);
% xw_ub(1:config.n_v) = ones(config.n_v,1);


% Initialization parameters and variables 
%debug
%xw_pri = [0.001; 0.09; 0.0905];
 %xw_pri_t = [0.0011; 0.0808; 0.0838]; % xw_pri = [0.0015; 0.0902; 0.0907];




Nk = 1;
 factor = [ones(1, Nk), zeros(1, config.N-Nk)]';
ym_init_mx = [ym_store(:,1), zeros(plant.ny,config.N-1)] ;
if tsMHE_estimator==2 || (tsMHE_estimator==0 && onoffN==1)
ym_init_mx = [ym_store(:,1:2), zeros(plant.ny,config.N-1)] ;
end
ym_init = reshape(ym_init_mx, [],1);



val_init = zeros(po.n_var,1);
if tsMHE_estimator==1  % ts MHE
z_pri = config.Txz(ym_store(1,1), xw_pri(2:end) );
par_init = [z_pri; ym_init; Nk; factor;u_const];  % t_0=1;
val_init(1:config.n_x) = reshape( kron(ones(config.N+1,1), z_pri ),[],1); 
val_init( config.n_x+1:config.n_x+config.n_d)= reshape(kron(ones(config.N,1), d_est(2:end) ), [], 1 );
val_init(config.n_x+config.n_d+1:end) = reshape( kron(ones(config.N,1), config.h_dyn(z_pri, d_est(2:end), xw_pri(1)+d_est(2) ) ), [], 1 );
elseif tsMHE_estimator ==2
par_init = [xw_pri; ym_init; Nk; factor;u_const];  % t_0=1;
val_init(1:config.n_x) = reshape( kron(ones(config.N+1,1), xw_pri),[],1); 
val_init( config.n_x+1:config.n_x+config.n_d)= reshape(kron(ones(config.N+1,1), d_est ), [], 1 ); 
val_init(config.n_x+config.n_d+1:end) = reshape(kron(ones(config.N+1,1), plant.y_fcn_full_MHE(xw_pri,d_est) ), [], 1);
else
  par_init = [xw_pri; ym_init; Nk; factor;u_const];  % t_0=1;
val_init(1:config.n_x) = reshape( kron(ones(config.N+1,1), xw_pri),[],1); 
val_init( config.n_x+1:config.n_x+config.n_d)= reshape(kron(ones(config.N+1*onoffN,1), d_est ), [], 1 ); 
val_init(config.n_x+config.n_d+1:end) = reshape(kron(ones(config.N+1*onoffN,1), plant.y_fcn_full_MHE(xw_pri,d_est) ), [], 1);  
end




% store results
sol_store = zeros(po.n_var,config.N_sim);
x_pri_store = zeros(config.nx,config.N_sim);
est_store = zeros(plant.n+plant.nw,config.N_sim);
est_store_z = zeros(config.nx, config.N_sim);
err = zeros(1, config.N_sim);
if tsMHE_estimator==1
    est_store_z(:,1) = z_pri;
    [est_store(:,1)] = xw_pri;   %proj_zx(est_store_z(:,1), ym_store(1,:), config.T_mx, config.xlim);
    x_pri_store(:,1) = z_pri;
else
    est_store(:,1) = xw_pri;
    x_pri_store(:,1) = est_store(:,1);
end

cost_store = zeros(1,config.N_sim);
par_store = zeros(length(par_init)-config.N,config.N_sim);




%% 


Tsim = zeros(1,config.N_sim);

Iter = zeros(1,config.N_sim);
%disp('***************Simulation Started!********************');
for k= 2: config.N_sim  % debug
 tic

 %disp('Simulation step ' + string(k) + '/' + string(config.N_sim));
    sol = solver('x0' , val_init,...         % solution guess
                 'p', par_init,...         % parameters
                 'lbx', xw_lb,...           % lower bound on x
                 'ubx', xw_ub,...           % upper bound on x
                 'lbg', c_lb,...           % lower bound on g
                 'ubg', c_ub);             % upper bound on g
 
    if ~solver.stats().success
        error('Solver did not find optimal solution at step k = ' + string(k));
    else

    
 
     val_opt = full(sol.x);
    % extract optimal state,disturbance 
    x_opt = reshape(val_opt(1:config.n_x),[config.nx,config.N+1]);

      Nk = min(k, config.N);
    if tsMHE_estimator==1  % ts-MHE
        est_store_z(:, k) = x_opt(:, min(k, config.N+1));
        xw_pri_temp = est_store_z(:, k-Nk+1);
        [est_store(:,k), err(:,k)]=proj_zx(est_store_z(:,k), ym_store(:,k), config.T_mx, config.xlim);  %%%%%%%
    else
     est_store(:,k) = x_opt(:, min(k, config.N+1));
    xw_pri_temp = est_store(:, k-Nk+1);

    end
    
   %         Tsim(k) = solver.stats().t_proc_total;
        Tsim(k) = toc;
        inf_cons = solver.stats().iterations.inf_pr; % primal feasibility
        inf_dual = solver.stats().iterations.inf_du;  % dual feasibility
%         Iter(k) = solver.stats().iter_count;
    end

    cost_store(1,k) = full(sol.f);
    sol_store(:,k) = val_opt;
    x_pri_store(:,k)= xw_pri_temp;
  
    if tsMHE_estimator==2 || (onoffN==1 && tsMHE_estimator==0)
        kt = min(k, config.N_sim-1);
    ym_init_mx = [ym_store(:,k+1-Nk: kt+1), zeros(plant.ny,config.N-Nk)] ;  
    else
    ym_init_mx = [ym_store(:,k+1-Nk: k), zeros(plant.ny,config.N-Nk)] ;
    end
    ym_init = reshape(ym_init_mx, [],1); 
    factor = [ones(1, Nk), zeros(1, config.N-Nk)]';


   par_init = [xw_pri_temp; ym_init; Nk; factor;u_const]; 

   % update guess of variable for next step (optional)
    val_init= val_opt;   % deactive for mEMPC_controller=0, otherwise, infeasibility arises!
   


end 


%% de-scaling
x_real = real_traj(:,1:config.N_sim);
x_est = est_store;
ym_real =ym_store;






%%

tMul= sum(Tsim);
% %iMul = sum(Iter);
%disp('Simulation Time ' + string(tMul) + 's');
% % disp('Iteration ' + string(iMul) );
res.x_init = xw_init;

res.sol = sol_store;
res.cost = cost_store;  
res.ym_store = ym_real;
res.traj = x_real;
res.traj_est = x_est;
res.config = config;
res.tMul =  tMul;
res.Design_Data = res_name;

