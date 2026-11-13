
config.Time = 0.5; %[h]
config.N=config.Time*60*60/plant.dt;  % [s]/[s]
config.N_sim = 10; %[h]
config.cost_w = 25;
config.cost = @(w,u, Q) (-Q*w'*w+u'*u );
config.wx_init = [0.1; 15; 5e-4 ];

options = struct;
options.ipopt.max_iter = 500;
options.ipopt.print_level = 0;
options.print_time = 0;

config.n_v = config.N*plant.nu;
config.n_x =(config.N+1)*plant.n;
config.n_w = (config.N+1)*plant.nw;


