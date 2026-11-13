%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for the numerical example in the paper "strong nonlinear detectability and moving horizon estimation for nonlinear sytems with unknown inputs"
% Author: Yang Guo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all;
addpath('..\..\..\casadi-3.6.3-windows64-matlab2018b');  % path of installed casadi toolbox
addpath('MHE\');
addpath('RCCM_syn\');
addpath('Data_mat\');
addpath('oolbox_share\'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type_MHE ={'benchmarkMHE', 'tsMHE', 'MHE_UIE'}; %tsMHE: UISE using residual model,  MHE_UIE: UISE using full model

ioss_onoff = 0;  % verifcation of detectability, 0:no, 1:yes
v_noise = 1e-2;   %bound of measurement noises 
sim_N = 1;      % number of simulation runs
onoffN=0;  %setting for benchmark MHE: 0=use past measurements for, 1=use current measurments 
%%% setting for plotting
printing = 1;

res_name = {'Data_mat\Data_design_full_model.mat'}; % -->paper
% res_name = {'Data_mat\Data_design_full_Henten_extend_v4_mod_test.mat'};

if v_noise==0
sim_res_name = {'Data_mat\Sim_result_N30_v0_test'};  %---> paper noise-free
  
else
sim_res_name = {'Data_mat\Sim_result_N30_v3e6_test'};  %---> paper  with noise
end

%% load plant

      load_model_full
       save(res_name{1}, "plant","-append"); %%

%% Verification of (partial) detectability
%% debug
if ioss_onoff
  
   rccm_s = ioss_synthesis_grid(rccm_s); 
   rccm = ioss_synthesis_grid(rccm);
   rccm_f = ioss_synthesis_grid_f(rccm_f);
   save(res_name{1}, "rccm","rccm_s");
   save(res_name{1}, "rccm_f","-append"); %%
 
end
   
%%
  
load(res_name{1});  % 3 input

%% MHE config (Casadi)
import casadi.*

config.N= 30;  % maxmial horizon which gives optimal solutions

config.Time_Sim=30; %[s]
config.N_sim = config.Time_Sim/plant.dt+1; %[h]-->steps 


sim_vec  = zeros(2,sim_N); 
sim_evec = zeros(2, sim_N);

if 1
%% simulate benchmark MHE
tsMHE_estimator=0;

load_MHE_config_scaling
% build MHE
if onoffN 
[solver,po] = setup_MHE_UIE( config, options); 
else
  [solver,po] = setup_MHE( config, options);   
end
for ii = 1:sim_N
run_simulation
  sim_vec(tsMHE_estimator+1,ii)=tMul/(config.N_sim-1);  % average time for pro sample in ii-th run
  e_vec = vecnorm(x_real-x_est);
  e_sse = norm(e_vec.^2); % sum of square errors
  sim_evec(tsMHE_estimator+1,ii) = e_sse; 
  disp(['simulation: ', num2str(ii),'/',num2str(3*sim_N),' is finished!']);
disp('-------------------------------------------------');
end
save([sim_res_name{1},type_MHE{1},'.mat'], "res");

%% simulate ts_MHE
tsMHE_estimator=1;
load_MHE_config_scaling
% build MHE
[solver,po] = setup_MHE( config, options); 
for ii = 1:sim_N
run_simulation
  sim_vec(tsMHE_estimator+1,ii)=tMul/(config.N_sim-1);
    e_vec = vecnorm(x_real-x_est);
  e_sse = norm(e_vec.^2); % sum of square errors
  sim_evec(tsMHE_estimator+1,ii) = e_sse; 
    disp(['simulation: ', num2str(ii+sim_N),'/',num2str(3*sim_N),' is finished!']);
disp('-------------------------------------------------');
end
save([sim_res_name{1},type_MHE{2},'.mat'], "res");

end
%% simulate MHE-based UIE
tsMHE_estimator=2;
load_MHE_config_scaling
% build MHE
[solver,po] = setup_MHE_UIE( config, options); 
for ii = 1:sim_N
    run_simulation
% run_simulation_UIE
  sim_vec(tsMHE_estimator+1,ii)=tMul/(config.N_sim-1);
    e_vec = vecnorm(x_real-x_est);
  e_sse = norm(e_vec.^2); % sum of square errors
  sim_evec(tsMHE_estimator+1,ii) = e_sse; 
    disp(['simulation: ', num2str(ii+2*sim_N),'/',num2str(3*sim_N),' is finished!']);
disp('-------------------------------------------------');
end
save([sim_res_name{1},type_MHE{3},'.mat'], "res");

%%% compute max simulation time for sim_N runs. 
 sim_max = max(sim_vec,[],2);
  sse_max_all = max(sim_evec,[],2);
%%save simulation results
sim_max_pro = sim_max(1);
sse_max = sse_max_all(1);
save([sim_res_name{1},type_MHE{1},'.mat'],'sim_max_pro', 'sse_max', 'sim_N', '-append'); 
sim_max_pro = sim_max(2);
sse_max = sse_max_all(2);
save([sim_res_name{1},type_MHE{2},'.mat'],'sim_max_pro', 'sse_max', 'sim_N', '-append'); 
sim_max_pro = sim_max(3);
sse_max = sse_max_all(3);
save([sim_res_name{1},type_MHE{3},'.mat'],'sim_max_pro', 'sse_max', 'sim_N', '-append'); 


 %% plot simulation results
 if 1
%%% Plot setup
lw = 2;         % line width
lw2 = 3;
fs_labels = 18; % fontsize labels
fs_axis = 22;   % fontsize axis
fs_legend = 18;
fs_shift_y = [0.02 0 0];
fs_shift_x = [0 0 0];


psize = [13.65, 9.6];
ppos  = [0.3, 0.2, psize(1)-0.3 psize(2)-0.2];
gcaPos = [0.15 0.16 0.8 0.8];
Ftsize_axis = 22;
FTsize_axisnumbers = 18;
FTsize_legend=10;


folder='Sim_plots/';
%%%%%%%%%%%%%%%
  
% load simulaton results
 load([sim_res_name{1},type_MHE{1},'.mat']);
 real_traj = res.traj;
 real_traj_est = res.traj_est;

 load([sim_res_name{1},type_MHE{2},'.mat']);
 real_traj_m = res.traj;
 real_traj_est_m = res.traj_est;


 load([sim_res_name{1},type_MHE{3},'.mat']);
 real_traj_f = res.traj;
 real_traj_est_f = res.traj_est;

 close all
 cm=f_getColorMap;


% plot results
time= 0: length(real_traj(1,:))-1;


figure(11)  
plot(time.*0.25, (real_traj(2,:)),'Color',cm(14,:), 'LineWidth', lw);
hold on
plot(time.*0.25, (real_traj_est(2,:)),'Color',cm(5,:), 'LineStyle','--','LineWidth', lw);
plot(time.*0.25, (real_traj_est_m(2,:)),'Color',cm(1,:), 'LineStyle',':','LineWidth', lw2);
plot(time.*0.25, (real_traj_est_f(2,:)),'Color',cm(3,:), 'LineStyle',':','LineWidth', lw2);
set(gcf,'PaperUnits','centimeters','PaperSize',psize,'PaperPosition',ppos);
set(gca,'Position',gcaPos,'fontsize',FTsize_axisnumbers,'TickLabelInterpreter','latex');
h23_2=xlabel('Time','Interpreter','latex','FontSize', fs_labels);
h23_1= ylabel('$x^b$', 'Interpreter','latex','FontSize', fs_labels); % 'Position', [Xlb_y Ylb_y]
set(h23_1, 'Position',h23_1.Position+fs_shift_y);
set(h23_2, 'Position',h23_2.Position+fs_shift_x);

if printing
      name = sprintf('X_d_1');
    print([folder name],'-depsc','-painters')
    print([folder name],'-dpdf','-painters')
    print([folder name],'-dbmp')
    savefig([folder name])  
end


% 
figure(33) 
plot(time.*0.25, (real_traj(1,:)),'Color',cm(14,:), 'LineWidth', lw);
hold on
plot(time.*0.25, (real_traj_est(1,:)),'Color',cm(5,:), 'LineStyle','--','LineWidth', lw);
plot(time.*0.25, (real_traj_est_m(1,:)),'Color',cm(1,:), 'LineStyle',':','LineWidth', lw2);
plot(time.*0.25, (real_traj_est_f(1,:)),'Color',cm(3,:), 'LineStyle',':','LineWidth', lw2);
set(gcf,'PaperUnits','centimeters','PaperSize',psize,'PaperPosition',ppos)
set(gca,'Position',gcaPos,'fontsize',FTsize_axisnumbers,'TickLabelInterpreter','latex');
ylabel('$x^a$', 'Interpreter','latex','FontSize', fs_labels);


figure(44)
plot(time.*0.25, res.config.w_true,'Color',cm(14,:), 'LineWidth', lw);
set(gcf,'PaperUnits','centimeters','PaperSize',psize,'PaperPosition',ppos)
set(gca,'Position',gcaPos,'fontsize',FTsize_axisnumbers,'TickLabelInterpreter','latex');
h33_4=xlabel('Time','Interpreter','latex','FontSize', fs_labels);
h23_4=ylabel('$w$', 'Interpreter','latex','FontSize', fs_labels);
set(h33_4, 'Position',h33_4.Position+fs_shift_x);
set(h23_4, 'Position',h23_4.Position+fs_shift_y);


if printing
      name = sprintf('W_d');
    print([folder name],'-depsc','-painters')
    print([folder name],'-dpdf','-painters')
    print([folder name],'-dbmp')
    savefig([folder name])  
end


psize = [13.65, 9.6];
ppos  = [0.3, 0.2, psize(1)-0.3 psize(2)-0.2];

gcaPos = [0.15 0.16 0.8 0.8];
Ftsize_axis = 22;
FTsize_axisnumbers = 18;
FTsize_legend=10;

%%
figure(55)
plot(time.*0.25, vecnorm((real_traj_est(2:end,:)-real_traj(2:end,:))),'Color',cm(5,:), 'LineStyle',':', 'LineWidth', lw2);
hold on
plot(time.*0.25, vecnorm((real_traj_est_m(2:end,:)-real_traj(2:end,:))),'Color',cm(1,:), 'LineStyle',':', 'LineWidth', lw2);
plot(time.*0.25,vecnorm((real_traj_est_f(2:end,:)-real_traj(2:end,:))),'Color',cm(3,:), 'LineStyle',':', 'LineWidth', lw2);
set(gcf,'PaperUnits','centimeters','PaperSize',psize,'PaperPosition',ppos)
set(gca,'Position',gcaPos,'fontsize',FTsize_axisnumbers,'TickLabelInterpreter','latex', 'YScale', 'log');


h33_5=xlabel('Time','Interpreter','latex','FontSize', fs_labels);
ylabel('$ \|x-\widehat{x}\|$', 'Interpreter','latex','FontSize', fs_labels);
ylim([10e-5,1]);
set(h33_5, 'Position',h33_5.Position+fs_shift_x);


if 1 %printing
    name = sprintf('Err_norm');
    print([folder name],'-depsc','-painters')
    print([folder name],'-dpdf','-painters')
    print([folder name],'-dbmp')
    savefig([folder name])  
end




 end


